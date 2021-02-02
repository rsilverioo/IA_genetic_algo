# IA - Genetic algorithms (LibGA)

In this project, the LibGA library developed by Arthur L. Corcoran has been implemented to find the max clique of some common instances. It is considered an introductory project to genetic algorithms and is not intended to be overly complex or precise.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

What things you need to install the software and how to install them:

* Nothing at all, everything you need is in this repository.

### Usage example
The main configuration for LibGA is *GAconfig* and the rest is in the *configs* folder, which is where the bash scripts will take them from to change them automatically. 

The available instances are in the *instances* folder, but you can add any other instance you want. Be sure to add them in the header of the config file for the program to use.

The code is divided into 3 different programmes:

**GA Test:** which is the normal program for using LibGA and uses the configuration of the *GAconfig* file. It can be used with its bash file to automate the change of configurations with those in the *config* folder and perform 10 tests per configuration file. A .csv file is created with all results and also the individual results of each simulation are saved in a .txt.

**GA Test Report:** modification of GA Test in which only the *GAconfig* configuration is used. The value of *mu_rate* is progressively changed, with constant *x_rate*, and after *x_rate*, keeping *mu_rate* constant. It is recommended to use the bash script with the same name to perform the 20 repetitions (50 samples of each) automatically. A .csv file is created with the final results.

**GA Total Test:** exactly the same process as in GA Test Report but a simulation is performed by varying *mu_rate* and *x_rate* at the same time. In such a way that a 20x20 matrix is obtained with the best solution for each pair of *mu_rate* and *x_rate* values, which is saved in a .csv file. Example below for the instance p_hat300_1.

![](https://i.ibb.co/Kbbx81F/matrix.png)

## Contributing

1. Fork it (<https://github.com/rsilverioo/IA_genetic_algo/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

## Authors

* **Arthur L. Corcoran** - LibGA creator - [artcorcoran](https://github.com/artcorcoran)
* **Rodrigo Silverio** - Code author - [rsilverioo](https://github.com/rsilverioo)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Claudio Rossi - For minor changes in LibGA lib

