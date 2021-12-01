## Pattern recognition tasks
Tasks 6 and 5 only.

### Task 6:

1. Model a random sample of n-dimension objects.

   Modelled sample objects are distributed normally, the mean and the variance have to be predefined in a configuration file.

2. Plot objects' 2d-projections.

   A number of each dimension is defined in an application form (starting with 0).

3. Classify objects using chosen algorithm (Parzen windows, k nearest neighbours, Bayes classifier, Parametric classifier).

4. Calculate a classification error.

5. Calculate a transformation matrix.

   An element t<sub>ij</sub> of this matrix stands for the number of elements of a class j classified as a class j. 

You'll need a `test.txt` to define the number of classes, the number of dimensions, prior probabilities, for each class its distribution parameters. Using QtCreator open `lab6.pro`, run the project. Input the size of a training sample in a textbox signed `N train`, the size of a test sample in a `N test` textbox, click `Model` to model the sample, click `test` to test a chosen classifier.

![img](/img/task6.jpg)



### Task 5 (Bayesian classifier):

1. Model a random sample of n-dimension objects.

   Modelled sample objects are distributed normally, the mean and the variance have to be predefined in a configuration file.

2. Plot objects' 2d-projections.

   A number of each dimension is defined in an application form (starting with 0).

3. Classify objects to minimize risk and to minimize error.

4. Calculate an average risk and a classification error.

5. Calculate a transformation matrix.

   An element t<sub>ij</sub> of this matrix stands for the number of elements of a class j classified as a class j. 

You'll need a `config.txt` to define the number of classes, the number of dimensions, loss matrix, prior probabilities, for each class its distribution parameters. Using QtCreator open `lab5.pro`, run the project. Input the size of a sample in a textbox signed `N`, click `Model`.

The first plot is min error classification, the second is modelled data, the third plot is min risk classification, colors represent classes.

![img](/img/task5.jpg)













