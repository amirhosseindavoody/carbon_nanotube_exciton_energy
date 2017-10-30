// #include <stdio.h>
// #include <iostream>
// #include <string>
// #include <fstream>
//
// void write_log(std::string input)
// {
//
// 	std::ofstream log_file;
// 	log_file.open("log.dat", std::ios::app);
// 	log_file << input << std::endl;
//
// 	if (log_file.fail())
// 	{
// 		std::cout << "error in writing to the log file!!!" << std::endl;
// 		exit(EXIT_FAILURE);
// 	}
//
// 	log_file.close();
//
// 	std::cout << input << std::endl;
// }
