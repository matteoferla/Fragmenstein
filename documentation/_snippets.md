These are snippets that need to go somewhere
# 1

To make Victor carp at any exception, do

    Victor.error_to_catch = ()
    
This is because `try`/`except` can accept multiple error classes as a tuple.

# 2

`ConnectionError` is meant for comms connection errors, but here it is raised as bad bonding.

# 3

Show the steps

    from IPython.display import display
    
    monster = Monster([fore, aft]).merge()
    for m in monster.modifications:
        display(m)