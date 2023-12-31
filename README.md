# Chem 274A: Final Project
## Part 1 - C++ Matrix Class
### C++ Makefile
- `run` - builds and executes the matrix class executable, demonstrating its capabilities
- `clean` - clears any pre-existing executable/output files

To begin designing a matrix class in C++, working incrementally was the name of the game. By looking over the specifications and grouping certain features and methods together, I was able to create a plan of action to properly implement the matrix class, beginning with private members.  

1. Private members  
The private members of my matrix class include:
    - `n_rows` - the number of rows in a matrix
    - `n_columns` - the number of columns in a matrix
    - `data_` - a vector of templated vectors containing the data of a matrix

    Implementing a matrix this way and encapsulating the data within a vector of vectors comes with many advantages. Vectors are dynamically allocated, meaning memory management and size allocation are done automatically, reducing the risk of memory leaks. Vectors also come with functions that come into play in later implementations, such as `.resize()` and `.clear()`. Additionally, using a vector of vectors allows for the access of elements via both the row and column number in tandem, rather than using a flat, 1D vector. This change makes accessing elements more intuitive overall. Additionally, the entire class is templated, meaning the desired data types of `int`, `double`, and `complex` are all compatible and able to be stored within a matrix.

2. Constructors  
The three implemented constructors are as follows:
    - Default constructor with no arguments
    - Constructor with desired number of rows and columns
    - Copy constructor for an existing matrix

    Implementing the constructors was relatively easy, only having to consider what the starting parameters were for each case. The default constructor was constructed as a 2x2 matrix, for better visibility when printed. The second constructor takes in arguments for the number of rows and columns desired and creates an empty matrix based on those values. Lastly, the copy constructor takes in a reference to a constant reference to another matrix in order to create the copy.

3. Accessing and resizing  
    To access individual elements, the parenthesees operator was overloaded to index to the desired row and column. Since the data is templated, the return type of the element is `T`. Both a const and non-const of this operator are programmed, in order to be compatible with const matrices. The non-const version allows for re-assignment of an individual element (i.e. `matrix(0,0) = 1`), whereas the const version is a read-only operation. Matrix resizing is done via the `resize` member function, which takes in parameters for the new rows and new columns. The private variables are updated to these new values and, prior to resizing the `data_` vector, all pre-existing data is cleared, as instructed by the project specifications. That way, no pre-existing data is able to interfere with operations on a resized vector. The `data_` vector is then resized according to the new dimensions, with default values of whatever type (typically 0). It is important to note that there is no return for resizing, as the resizing operation is done in-place on the current matrix.

4. Mathematical matrix operations  
The implemented mathematical operations are as follows:
    - Matrix addition (via the `+` operator)
    - Matrix subtraction (via the `-` operator)
    - Matrix multiplication (via the `*` operator)
    - Matrix element-wise multiplication

    Development of these operations required some mental refreshers on topics related to linear algebra. Matrix addition, subtraction, and element-wise multiplication all required the dimensions of both the left hand and right hand side matrices to have the same dimensions. To enforce these conditions, an `if... throw` statement is included in each function that throws an `invalid_argument` error if the rows or columns are unequal. Similarly, the matrix multiplication has its own enforced condition, being the number of columns in the left matrix must be equal to the number of rows in the right matrix. These operations are const correct and return a new matrix as a result, which is constructed by iterating through each element position and assigning it a value based on each mathematical operation.

5. Matrix re-assignment and filling  
    Matrix re-assignment is done via the `=` operator, which has been overloaded to take in either another matrix or a scalar. To re-assign a matrix to another matrix, it is first resized to the dimensions of the right hand matrix, then iterates through each position and updates the element. When assigning a matrix to a scalar value, the operation uses the overloaded `=` operator function to fill all element positions with that scalar value. A similar procedure is done using the `.fill()` method. All these methods are done in-place and return the current matrix by reference.

6. Matrix comparison  
    For matrix comparison, the `==` operator is overloaded to return a boolean value. In order to compare two matrices together, their dimensions are first checked prior to iterating. If they are equal, then `std::equal` is utilized to iterate through both matrices and their elements. Using a standard library function prevents having to iterate through each element and comparing the data.

7. Matrix printing  
    To print a matrix to an output stream object, some research was done on how to work with the `friend` declaration. This keyword gives a function access to another class's private member variables, which is required to feed each element into the output stream via `m.data[i][j]`. A singular whitespace follows each element in a row, and a singular line skip follows at the end of a row, leading to the clean printing of a matrix. By overloading the output stream operator, `<<`, the matrix can be accessed via `std::cout` and be printed. An extra print function is also provided that reuses the `<<` operator.

8. Linear algebra and complex number operations  
The last features implemented in the matrix class were to return:
    - The transpose of a matrix
    - The complex conjugate of a matrix
    - The conjugate transpose of a matrix
    - In-place versions of each of the above operations

    For the non in-place versions of these functions, a new matrix is returned with each element being changed accordingly. In the transpose function, a new matrix is constructed with opposite columns and rows, with each element being assigned to its transpose. The complex conjugate of a matrix utilizes the `std::conj` function from the `complex` library to return the complex conjugate of an element. The conjugate transpose takes advantage of the pre-existing functions and assigns a temporary matrix that first returns the complex conjugate, then assigns the returned matrix to the transpose of the temporary matrix. In doing so, the code is re-used rather than copied, leading to cleaner code. The in-place versions of each function do not return anything, but instead re-assign the current matrix to the corresponding function via `*this`. Again, code is re-used to save time and promote efficiency.

9. Debugging  
    Valgrind and address sanitizers were used in tandem to check for any issues in memory. Thankfully, due to the use of vectors, memory allocation was mostly an afterthought and was handled automatically by the compiler. Common bugs that appeared during development were often due to const correctness and making sure that certain functions were const correct (i.e. mathematical operations) and others weren't (i.e. element re-assignment, in-place operations).

10. Tests  
    The written test file demonstrates the capabilities of the class, including all of the above operations. The output contains many examples of matrices with different datatypes, the in-place and non in-place operations, and the printing of these matrices. Printing to an output file is also included to show that the file stream is compatible with the overloaded `<<` operator.

Overall, this assignment was a great exercise into class design and what goes into implementing a container that is not native to the language. Const correctness, operator overloading, error handling, and in-place implementation were major players in making the class as robust as it is. Some further improvements could be to include functions to calculate eigenvalues/eigenvectors of a given matrix, as well as other linear algebra operations. Regardless, this assignment encapsulated the lessons taught over the course of the semester and rewarded me with a better, fundamental understanding of C++ and its nuances.



## Part 2 - Python Molecule Class and Substructure Searching
### Python Makefile
- `environment` - builds the environment with all required packages/libraries for the molecule class
- `test` - runs pytest test cases located in `test_mol.py`

To begin implementing a molecule class, capable of substructure searching, it was necessary to learn what a molecular fingerprint is and what it entails. Essentially, a molecular fingerprint is a vector of bits that represents the functional groups present and/or absent in a certain molecule. While seemingly abstract and unclear at first, proper planning and trial-and-error came into play when being able to properly implement this class and its capabilities.

1. Initialization and drawing  
    One of the first noticeable things about this assignment was that NetworkX was encouraged for the development of a molecule class. Additionally, the `parse_sdf()` function is provided, reminiscent of Problem Set 3. Initially, I had to distinguish the difference between composition and inheritance, prior to considering them in development. *Composition* is "has a" relationship between classes, meaning a molecule class would have a relationship with the `MolGraph` class from Problem Set 3. *Inheritance* is a "is a" relationship between classes, meaning a molecule class would be the child class of the `MolGraph` class. After toying around with initial implementations, I chose to rewrite the class entirely, but utilize the concepts and members the original `MolGraph` class contained. This is because the provided `parse_sdf()` function works a bit differently, requiring a new implementation altogether. The class still takes in a path to a `.sdf` file, parses through that file, and creates a NetworkX graph according to the file's contents. However, edge_labels are added to better visualize the bond order between two atoms, which is necessary for molecular fingerprinting. The `draw_graph()` function is mostly the same as before, but does not include a color map, as the previous color map only included few elements (HCNO). Given that the main role of this class is for substructure searching, I figured it would be sufficient to have the graphical representation of a molecule be more barebones, since we are more interested in the bond paths and bond orders between atoms.

2. Special method: `_ipython_display_`  
    One requirement of this class was to provide a graphical representation of a molecule specifically in IPython/Jupyter Notebooks. In these notebooks, a special method, `_ipython_display_` can be overloaded to print a visual if the last line of a code cell is an object. By calling `draw_graph()` in this method, we are able to see a visual of a molecule more easily in these notebooks, which is helpful in identifying potential functional groups.

3. Molecular fingerprinting  
    The bulk of this assignment was figuring out how to implement the molecular fingerprint of a molecule. This was a new topic for me and you can only take away so much from 2 pages of specifications. Thankfully, pseudocode was provided that helped simplify the implementation of the actual fingerprint, but what about the paths? This was the challenge for me, and I headed into the NetworkX documentation. After searching around, I arrived at the function `all_simple_paths`, which generated all paths in a graph starting at a node and ending at a target node, without repeating nodes and up to a cutoff depth. This was the exact function I needed since it did not repeat nodes and could be customized to end at a certain depth of 6, unlike the implementation in RDKit. Each path was processed to its corresponding SMILES equivalent, which could then be hashed. A fingerprint was interpreted as a NumPy array of 1024 zeroes (compared to RDKit's 2048 bit vector), which has two random positions set to 1 per hashed path. This implementation was a great mental exercise in trial-and-error, debugging, and creativity. The corresponding `fingerprint()` method was decorated with a `property` decorator, allowing it to be accessed as an object member for subsequent methods.

4. Equivalence overloading  
    To access and overload the equivalence operator, the dunder method, `__eq__`, was modified to first check the right hand side for its type via `isinstance`. If the `other` parameter was not a molecule, then a `TypeError` would be raised. The equivalence operator returns whether or not the fingerprints of both molecules are equal to each other, aka if all the positions in the array are the same. This seemed to be the best option for verifying equivalence as it utilized pre-existing code and features rather than comparing contents from `.sdf` files.

5. Substructure searching  
    To check a substructure's presence in a given molecule, I first did the same type-checking procedure done in the `__eq__` method, with a molecule object in the `.sdf` version versus a string in the SMILES version. `check_substructure_sdf` takes in another molecule object and calculates the difference between the two fingerprints. Clever implementation of the fingerprint as a NumPy array allowed for the use of NumPy array operations, such as array subtraction. Given that random positions in the array would be assigned to a value of 1, a matching substructure would return a matching array with all values being greater than or equal to 0, since the substructure fingerprint would only contain 1's in corresponding indices to the original fingerprint. Returning `np.all(match >=0)` verifies whether or not the resulting array has any negative values, meaning a given substructure is not found in a molecule. The same is done in the SMILES version, but the input is given as a SMILES path string, which is then hashed to create a fingerprint. This implementation, while not as tedious as fingerprinting, was another exercise in intuitive and creative thinking.

6. Testing  
    Initial tests and debugging of the source file were done in a developmental Jupyter Notebook, which was very helpful in seeing the outputs of each incremental implementation. Test cases were added to its own source file and run via pytest, which allowed for the parametrization of certain test cases, such as type-checking an instance and verifying whether a molecule has an equivalent or a substructure.

7. Sample Jupyter Notebook  
    The sample Jupyter Notebook, `mol_search.ipynb`, demonstrates the capabilities of the class and its substructure searching methods, using `.sdf` files. I figured this implementation would be overall more intuitive compared to command line searching, since a substructure screen in the command line would only be compatible with string inputs via `sys.argv` or `sys.argparse`. By implementing both methods for use in a Jupyter Notebook, a molecule can first be visualized properly via matplotlib and screened for certain functional groups, to ensure robustness of code. Overall, my class design was more tailored towards use in interactive notebooks since the Problem Set 3's class was very similar in concept. However, another approach could have been more command line based and oriented, in order to more quickly search for substructures.

Overall, the development of this class was much more challenging conceptually compared to the matrix class in C++, but more intellectually rewarding. It was another great exercise in class design, particularly in making decisions that ultimately the class's use(notebooks vs. command line). While composition neither inheritance were used in this project, both could have been certainly viable from Problem Set 3's class. Python still ends up being my preferred language for programming and development, but both parts of this final project have given me a greater sense of appreciation and understanding for both languages, which will no doubt be carried onwards in future classes.

Thank you Sifei, Mayank, Prof. Nash, and Prof. Pritchard for your help/work this semester!