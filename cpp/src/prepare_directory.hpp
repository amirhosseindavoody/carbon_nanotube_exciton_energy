#ifndef _prapare_directory_hpp_
#define _prapare_directory_hpp_

#include <experimental/filesystem>
#include <stdexcept>

inline std::experimental::filesystem::directory_entry prepare_directory(const std::string path, const bool keep_old_files=true)
{
  std::cout << "\n..." << std::endl;

  namespace fs = std::experimental::filesystem;
  fs::directory_entry directory;

  directory.assign(path);
  std::cout << "preparing directory: " << directory.path() << std::endl;

  if (not fs::exists(directory.path()))
  {
    // std::cout << "warning: output directory does NOT exist!!!" << std::endl;
    std::cout << "created the directory!" << std::endl;
    fs::create_directories(directory.path());
    if (not fs::is_directory(directory.path())) throw std::invalid_argument("The input value for output directory is not acceptable.");
    return directory;
  }

  if (fs::is_directory(directory.path()))
  {
    if (not fs::is_empty(directory.path()))
    {
      std::cout << "warning: output directory is NOT empty!!!" << std::endl;
      if (keep_old_files)
      {
        int count = 1;
        while (fs::exists(directory.path().string()+"."+std::to_string(count))){
          count++;
        }
        fs::path new_path = directory.path().string()+"."+std::to_string(count);
        std::cout << "renaming the existing directory to " << new_path << std::endl;
        fs::rename(directory.path(),new_path);
      } 
      else
      {
        std::cout << "deleting the existing directory!!!" << std::endl;
        fs::remove_all(directory.path());
      }
      fs::create_directories(directory.path());
    }
  }
  else
  {
    throw std::invalid_argument("The input value for output directory is not acceptable.");
  }

  std::cout << "...\n" << std::endl;

  return directory;

};

#endif //_prapare_directory_hpp_