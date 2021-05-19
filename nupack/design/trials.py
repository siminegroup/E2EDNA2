import random, datetime, warnings, sys

from tempfile import NamedTemporaryFile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from .components import TimeInterval, WriteToFileCheckpoint, StopCondition
import logging

# Overwrite this to silence or customize output
# e.g. log = lambda s: None

logger = logging.getLogger(__name__)
log = lambda s: logger.warning(s)

################################################################################

class Trial:
    '''Class representing a single trial of a design optimization'''

    def __init__(self, design, *, index, checkpoint, interval):
        self.design = design
        self.spec = design.spec.copy()
        self.spec.parameters.seed = random.randint(1,100000)
        self.future = None
        self.index = int(index)
        self.checkpoint = None if checkpoint is None else checkpoint/str(index)
        self.interval = interval
        self.condition = None
        self.output = None

    def setup_files(self):
        '''Set up filesystem storage'''
        self.spec.parameters.log = NamedTemporaryFile(delete=False, suffix='.txt').name

    def run(self, restart):
        '''Main task within the future: runs the design'''
        if restart is None:
            file, res = None, None
        elif isinstance(restart, (Path, str)):
            file, res = self._load_checkpoint(Path(restart)/'checkpoint')
        else:
            file, res = None, restart

        if self.checkpoint is None:
            chk = ' with no checkpoints'
        else:
            chk = ' with checkpoints in %s' % self.checkpoint

        if res is None:
            log('Starting design {} from scratch'.format(self.index) + chk)
        elif file is None:
            log('Resuming design {} from loaded result'.format(self.index) + chk)
        else:
            log('Resuming design {} from file {}'.format(self.index, file) + chk)

        if self.checkpoint is None:
            self.condition = StopCondition()
            result = self.spec(restart=None if res is None else res.raw,
                checkpoint_condition=self.condition)
        else:
            checkpoint = self.checkpoint/'checkpoint'
            self.checkpoint.mkdir(exist_ok=True, parents=True)

            self.condition = TimeInterval(self.interval)
            result = self.spec(restart=None if res is None else res.raw,
                checkpoint_condition=self.condition,
                checkpoint_handler=WriteToFileCheckpoint(checkpoint))

            self._final().write_text(result.to_json().dump(indent=4))

        result = self.design.build_result(result)
        self.output = result
        return result

    def check(self):
        '''Check that the future is present and not failed (non-blocking)'''
        if self.future is None:
            raise RuntimeError('design optimization has not started yet')
        if self.future.done() and self.future.exception() is not None:
            self.future.result() # throw current exception

    def cancel(self):
        '''Cancel the design if it is running'''
        if self.condition is not None:
            self.condition.cancel()
        if self.future is not None:
            self.future.cancel()

    def _final(self):
        return self.checkpoint/'final.json'

    def _checkpoints(self, where=None):
        if where is None:
            where = self.checkpoint/'checkpoint'
        if not where.exists():
            return []
        return sorted((i for i in where.iterdir() if i.suffix == '.json'),
            reverse=True, key=lambda p: datetime.datetime.fromisoformat(p.stem))

    def current_result(self):
        '''
        Load latest result from save directory
        Return file path and result object
        Return None for each of those if unavailable
        '''
        if self.future is None:
            return None, None

        self.check()

        if self.checkpoint is None:
            return None, self.output

        final = self._final()
        if final.exists():
            try:
                return final, self.design.load_result(final)
            except Exception as e:
                log('Failed to load final result: {}'.format(e))

        if not self.checkpoint.exists():
            return None, None

        return self._load_checkpoint()

    def final_result(self):
        if self.future is None:
            return None

        self.check()
        return self.output

    def _load_checkpoint(self, where=None):
        for path in self._checkpoints(where):
            try:
                res = self.design.load_result(path)
                return path, res
            except Exception as e:
                log('Failed to load checkpoint: {}'.format(e))
        return None, None

    def cleanup(self, keep=5):
        '''Remove all but the most recent `keep` checkpoint files'''
        if self.checkpoint is None or not self.checkpoint.exists():
            return
        for p in self._checkpoints()[keep:]:
            p.unlink()

################################################################################

class Optimization:
    '''
    Class representing a set of currently running design trials
    '''

    trial_class = Trial

    def __init__(self, design, trials=4, *, checkpoint=None, interval=600):
#         if checkpoint is None:
#             warnings.warn('''checkpoint is set to None.
# * No check points will be saved and intermediate results are disabled!
# * Use Optimization(design, checkpoint = \'some-directory/\') to save check points and enable intermediate results.''')

        self.started = False
        self.checkpoint = None if checkpoint is None else Path(checkpoint)
        self.trials = [self.trial_class(design, index=i, checkpoint=self.checkpoint, interval=interval) for i in range(trials)]
        self.pool = ThreadPoolExecutor(trials)

    def start(self, restart=None):
        '''Start the design running.'''
        if self.started:
            raise Exception('Designs already running')
        self.started = True

        if self.checkpoint is not None:
            self.checkpoint.mkdir(exist_ok=True, parents=True)

        if restart is None:
            restart = [None] * len(self.trials)
        elif isinstance(restart, (Path, str)):
            restart = [Path(restart)/str(i) for i in range(len(self.trials))]
        else:
            if len(restart) != len(self.trials):
                raise ValueError('Wrong number of restarts given')

        # if not all(r is None or isinstance(r, Result) for r in restart):
        #     raise TypeError('Expected instance of Result (or None)')

        for t, r in zip(self.trials, restart):
            t.cancel()
            t.setup_files()
            t.future = self.pool.submit(t.run, r)

    def cleanup(self, n=5):
        '''Deletes old output files (keeps only the n newest)'''
        for t in self.trials:
            t.cleanup(n)

    def stop(self, force=True):
        '''Stop the design. Does not delete any output files.'''
        if not self.started:
            raise RuntimeError('No designs are running')

        if self.checkpoint is None and not force:
            raise RuntimeError('''No checkpoint directory was specified.
No outputs will be generated for designs that have not completed.
If you are sure you want to terminate designs, use stop(force = True).''')

        self.started = False
        for t in self.trials:
            t.cancel()

    def wait(self):
        '''
        Load the final results of each design. Will block until each design is complete.
        '''
        return [t.future.result() for t in self.trials]

    def current_results(self):
        '''
        Return the current checkpoint result of every trial.
        '''
        return [t.current_result()[1] for t in self.trials]


    def final_results(self):
        '''
        Return the final checkpoint result or None for every trial.
        '''
        return [t.final_result() for t in self.trials]

    def __del__(self):
        if self.started:
            self.stop(force=True)

