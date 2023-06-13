## Classes calling classes

> :label: Note to avoid confusion. 
> Herein the terminology is being used strictly: 
> the term `class` is used to refer to a `class` proper and not an `instance` of that class.

The `Laboratory` class has a class attribute `.Victor`,
which is simply the class `Victor`, unless a subclass of `Victor` is passed.
This class will be called to create a `Victor` instance for each task within the `Laboratory` instance.
In `fragmenstein.faux_victors` are a few such subclasses of `Victor`, which were used for testing.

Likewise, `Victor` class has a class attribute `.Monster`, which is simply the class `Monster`.
This is not to be confunded with the `Monster` instance, `monster` that is created by `Victor` instantiation.

Until [issue 34 is resolved](https://github.com/matteoferla/Fragmenstein/issues/34),
there is a confusing detail in that `Victor.monster_throw_on_discard` is a class attribute 
that when `Monster` is instantiated sets the `Monster` instance's attribute `throw_on_discard`.

