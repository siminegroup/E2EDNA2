def example_function(home, hight=20, **kwargs):
    print(home)
    for key, value in kwargs.items():
        print(f"{key}: {value}")

# Call the function with undesignated keyword arguments
example_function('earth',name="John", age=25, city="New York")
