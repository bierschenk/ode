=======
History
=======


0.4.0 (2019-07-21)
------------------

* Now requires numpy
* Replaces list comprehension addition/multiplication with array math
* update docs/double-pendulum-example to use an array dot-function and
  place it on the usage-->examples page of the documentation.
* The date on 0.3.0 was wrong -- should have been 2019-07-20

0.3.0 (2019-06-10)
------------------

* BREAKS COMPATIBILITY
* Output changes from [t, [[x0, v0], [x1, v1] ...]]
  to [t, [[x0, x1, ...], [v0, v1, ...]]].

0.2.0 (2018-06-10)
------------------

* Refactor to improve readability and code reuse.
* Clean tempating and testing
* Update documentation


0.1.2 (2018-03-06)
------------------

* Add verlet integration


0.1.1 (2018-03-05)
------------------

* Add generator versions of Euler and backward Euler methods.


0.1.0 (2018-03-03)
------------------

* First release on PyPI.
