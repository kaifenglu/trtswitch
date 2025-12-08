#ifndef __DATAFRAME_LIST__
#define __DATAFRAME_LIST__

#include <iostream>
#include <cmath>
#include <numeric>
#include <vector>
#include <string>
#include <limits>
#include <map>
#include <unordered_map>
#include <variant>
#include <algorithm>
#include <functional>
#include <stdexcept>

#include <Rcpp.h>


// Forward declaration for the recursive type
struct ListCpp;


// A C++ DataFrame structure to manipulate data frames efficiently
struct DataFrameCpp {
  std::vector<std::string> names_;
  std::unordered_map<std::string, std::vector<double>> numeric_cols;
  std::unordered_map<std::string, std::vector<int>> int_cols;
  std::unordered_map<std::string, std::vector<bool>> bool_cols;
  std::unordered_map<std::string, std::vector<std::string>> string_cols;
  
  size_t nrows() const {
    if (!names_.empty()) {
      const std::string& name = names_.front();
      if (numeric_cols.count(name)) return numeric_cols.at(name).size();
      if (int_cols.count(name)) return int_cols.at(name).size();
      if (bool_cols.count(name)) return bool_cols.at(name).size();
      if (string_cols.count(name)) return string_cols.at(name).size();
    }
    return 0;
  }
  
  size_t size() const {
    return names_.size();
  }
  
  const std::vector<std::string>& names() const {
    return names_;
  }
  
  bool containElementNamed(const std::string& name) const {
    return numeric_cols.count(name) || int_cols.count(name) || 
      bool_cols.count(name) || string_cols.count(name);
  }
  
  void check_row_size(size_t size, const std::string& name) const {
    size_t current_rows = nrows();
    if (current_rows > 0 && size != current_rows)
      throw std::runtime_error("Column '" + name + "' has inconsistent number of rows");
  }
  
  void push_back(const std::vector<double>& col, const std::string& name) {
    if (containElementNamed(name)) 
      throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    numeric_cols[name] = col;
    names_.push_back(name);
  }
  
  void push_back(const std::vector<int>& col, const std::string& name) {
    if (containElementNamed(name))
      throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    int_cols[name] = col;
    names_.push_back(name);
  }
  
  void push_back(const std::vector<bool>& col, const std::string& name) {
    if (containElementNamed(name))
      throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    bool_cols[name] = col;
    names_.push_back(name);
  }
  
  void push_back(const std::vector<std::string>& col, const std::string& name) {
    if (containElementNamed(name))
      throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    string_cols[name] = col;
    names_.push_back(name);
  }
  
  /**
   * Pushes a double scalar to the back of the DataFrame by expanding it to a vector 
   * matching the DataFrame's current row count.
   */
  void push_back(double value, const std::string& name) {
    size_t current_rows = nrows();
    if (current_rows == 0 && !names_.empty()) {
      throw std::runtime_error("Cannot push scalar '" + name + "': DataFrame is currently empty (0 rows).");
    } else if (current_rows == 0) {
      // If DataFrame is empty, we can initialize it with one row
      current_rows = 1;
    }
    std::vector<double> expanded_col(current_rows, value);
    push_back(expanded_col, name); // Delegate to the vector overload
  }
  
  /**
   * Pushes an int scalar to the back of the DataFrame by expanding it to a vector 
   * matching the DataFrame's current row count.
   */
  void push_back(int value, const std::string& name) {
    size_t current_rows = nrows();
    if (current_rows == 0 && !names_.empty()) {
      throw std::runtime_error("Cannot push scalar '" + name + "': DataFrame is currently empty (0 rows).");
    } else if (current_rows == 0) {
      // If DataFrame is empty, we can initialize it with one row
      current_rows = 1;
    }
    std::vector<int> expanded_col(current_rows, value);
    push_back(expanded_col, name); // Delegate to the vector overload
  }
  
  /**
   * Pushes a bool scalar to the back of the DataFrame by expanding it to a vector 
   * matching the DataFrame's current row count.
   */
  void push_back(bool value, const std::string& name) {
    size_t current_rows = nrows();
    if (current_rows == 0 && !names_.empty()) {
      throw std::runtime_error("Cannot push scalar '" + name + "': DataFrame is currently empty (0 rows).");
    } else if (current_rows == 0) {
      // If DataFrame is empty, we can initialize it with one row
      current_rows = 1;
    }
    std::vector<bool> expanded_col(current_rows, value);
    push_back(expanded_col, name); // Delegate to the vector overload
  }
  
  /**
   * Pushes a string scalar to the back of the DataFrame by expanding it to a vector 
   * matching the DataFrame's current row count.
   */
  void push_back(const std::string& value, const std::string& name) {
    size_t current_rows = nrows();
    if (current_rows == 0 && !names_.empty()) {
      throw std::runtime_error("Cannot push scalar '" + name + "': DataFrame is currently empty (0 rows).");
    } else if (current_rows == 0) {
      // If DataFrame is empty, we can initialize it with one row
      current_rows = 1;
    }
    std::vector<std::string> expanded_col(current_rows, value);
    push_back(expanded_col, name); // Delegate to the vector overload
  }
  
  /**
   * Pushes a matrix (row-major std::vector<std::vector<double>>) to the back of the DataFrame.
   * Creates 1 column named 'base_name' if p==1, otherwise p columns named 'base_name.'1 to 'base_name.'p.
   */
  void push_back(const std::vector<std::vector<double>>& matrix_data_rows, const std::string& base_name) {
    size_t num_rows = matrix_data_rows.size(); // Number of rows in the matrix
    if (num_rows == 0) return; 
    
    size_t p = matrix_data_rows[0].size(); // Number of columns in the matrix (p)
    
    // Check consistency of column counts across rows
    for (size_t i = 1; i < num_rows; ++i) {
      if (matrix_data_rows[i].size() != p) {
        throw std::runtime_error("Input matrix has inconsistent column counts across rows.");
      }
    }
    
    // Check row size compatibility with existing DataFrame rows
    check_row_size(num_rows, base_name);
    
    if (p == 1) {
      // Case 1: p == 1, create one column named 'base_name' (e.g., "v")
      if (containElementNamed(base_name)) 
        throw std::runtime_error("Column '" + base_name + "' already exists.");
      
      // We need to extract the single value from each row into a flat vector
      std::vector<double> single_column(num_rows);
      for(size_t row = 0; row < num_rows; ++row) {
        single_column[row] = matrix_data_rows[row][0];
      }
      
      numeric_cols[base_name] = std::move(single_column);
      names_.push_back(base_name);
      
    } else {
      // Case 2: p > 1, create p columns named 'base_name'1 to 'base_name'p (e.g., "v1" to "vp")
      
      // First, we must transpose/reorganize the data from row-major to column-major for storage
      std::vector<std::vector<double>> columns_data(p, std::vector<double>(num_rows));
      
      for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < p; ++col) {
          columns_data[col][row] = matrix_data_rows[row][col];
        }
      }
      
      // Now push back the newly formatted columns
      for (size_t col = 0; col < p; ++col) {
        std::string col_name = base_name + "." + std::to_string(col + 1);
        
        if (containElementNamed(col_name)) 
          throw std::runtime_error("Column '" + col_name + "' already exists.");
        
        numeric_cols[col_name] = columns_data[col];
        names_.push_back(col_name);
      }
    }
  }
  
  /**
   * Pushes a matrix (row-major std::vector<std::vector<int>>) to the back of the DataFrame.
   * Creates 1 column named 'base_name' if p==1, otherwise p columns named 'base_name.'1 to 'base_name.'p.
   */
  void push_back(const std::vector<std::vector<int>>& matrix_data_rows, const std::string& base_name) {
    size_t num_rows = matrix_data_rows.size(); // Number of rows in the matrix
    if (num_rows == 0) return; 
    
    size_t p = matrix_data_rows[0].size(); // Number of columns in the matrix (p)
    
    // Check consistency of column counts across rows
    for (size_t i = 1; i < num_rows; ++i) {
      if (matrix_data_rows[i].size() != p) {
        throw std::runtime_error("Input matrix has inconsistent column counts across rows.");
      }
    }
    
    // Check row size compatibility with existing DataFrame rows
    check_row_size(num_rows, base_name);
    
    if (p == 1) {
      // Case 1: p == 1, create one column named 'base_name' (e.g., "v")
      if (containElementNamed(base_name)) 
        throw std::runtime_error("Column '" + base_name + "' already exists.");
      
      // We need to extract the single value from each row into a flat vector
      std::vector<int> single_column(num_rows);
      for(size_t row = 0; row < num_rows; ++row) {
        single_column[row] = matrix_data_rows[row][0];
      }
      
      int_cols[base_name] = std::move(single_column);
      names_.push_back(base_name);
      
    } else {
      // Case 2: p > 1, create p columns named 'base_name'1 to 'base_name'p (e.g., "v1" to "vp")
      
      // First, we must transpose/reorganize the data from row-major to column-major for storage
      std::vector<std::vector<int>> columns_data(p, std::vector<int>(num_rows));
      
      for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < p; ++col) {
          columns_data[col][row] = matrix_data_rows[row][col];
        }
      }
      
      // Now push back the newly formatted columns
      for (size_t col = 0; col < p; ++col) {
        std::string col_name = base_name + "." + std::to_string(col + 1);
        
        if (containElementNamed(col_name)) 
          throw std::runtime_error("Column '" + col_name + "' already exists.");
        
        int_cols[col_name] = columns_data[col];
        names_.push_back(col_name);
      }
    }
  }
  
  /**
   * Pushes a matrix (row-major std::vector<std::vector<bool>>) to the back of the DataFrame.
   * Creates 1 column named 'base_name' if p==1, otherwise p columns named 'base_name.'1 to 'base_name.'p.
   */
  void push_back(const std::vector<std::vector<bool>>& matrix_data_rows, const std::string& base_name) {
    size_t num_rows = matrix_data_rows.size(); // Number of rows in the matrix
    if (num_rows == 0) return; 
    
    size_t p = matrix_data_rows[0].size(); // Number of columns in the matrix (p)
    
    // Check consistency of column counts across rows
    for (size_t i = 1; i < num_rows; ++i) {
      if (matrix_data_rows[i].size() != p) {
        throw std::runtime_error("Input matrix has inconsistent column counts across rows.");
      }
    }
    
    // Check row size compatibility with existing DataFrame rows
    check_row_size(num_rows, base_name);
    
    if (p == 1) {
      // Case 1: p == 1, create one column named 'base_name' (e.g., "v")
      if (containElementNamed(base_name)) 
        throw std::runtime_error("Column '" + base_name + "' already exists.");
      
      // We need to extract the single value from each row into a flat vector
      std::vector<bool> single_column(num_rows);
      for(size_t row = 0; row < num_rows; ++row) {
        single_column[row] = matrix_data_rows[row][0];
      }
      
      bool_cols[base_name] = std::move(single_column);
      names_.push_back(base_name);
      
    } else {
      // Case 2: p > 1, create p columns named 'base_name'1 to 'base_name'p (e.g., "v1" to "vp")
      
      // First, we must transpose/reorganize the data from row-major to column-major for storage
      std::vector<std::vector<bool>> columns_data(p, std::vector<bool>(num_rows));
      
      for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < p; ++col) {
          columns_data[col][row] = matrix_data_rows[row][col];
        }
      }
      
      // Now push back the newly formatted columns
      for (size_t col = 0; col < p; ++col) {
        std::string col_name = base_name + "." + std::to_string(col + 1);
        
        if (containElementNamed(col_name)) 
          throw std::runtime_error("Column '" + col_name + "' already exists.");
        
        bool_cols[col_name] = columns_data[col];
        names_.push_back(col_name);
      }
    }
  }
  
  /**
   * Pushes a matrix (row-major std::vector<std::vector<std::string>>) to the back of the DataFrame.
   * Creates 1 column named 'base_name' if p==1, otherwise p columns named 'base_name.'1 to 'base_name.'p.
   */
  void push_back(const std::vector<std::vector<std::string>>& matrix_data_rows, const std::string& base_name) {
    size_t num_rows = matrix_data_rows.size(); // Number of rows in the matrix
    if (num_rows == 0) return; 
    
    size_t p = matrix_data_rows[0].size(); // Number of columns in the matrix (p)
    
    // Check consistency of column counts across rows
    for (size_t i = 1; i < num_rows; ++i) {
      if (matrix_data_rows[i].size() != p) {
        throw std::runtime_error("Input matrix has inconsistent column counts across rows.");
      }
    }
    
    // Check row size compatibility with existing DataFrame rows
    check_row_size(num_rows, base_name);
    
    if (p == 1) {
      // Case 1: p == 1, create one column named 'base_name' (e.g., "v")
      if (containElementNamed(base_name)) 
        throw std::runtime_error("Column '" + base_name + "' already exists.");
      
      // We need to extract the single value from each row into a flat vector
      std::vector<std::string> single_column(num_rows);
      for(size_t row = 0; row < num_rows; ++row) {
        single_column[row] = matrix_data_rows[row][0];
      }
      
      string_cols[base_name] = std::move(single_column);
      names_.push_back(base_name);
      
    } else {
      // Case 2: p > 1, create p columns named 'base_name'1 to 'base_name'p (e.g., "v1" to "vp")
      
      // First, we must transpose/reorganize the data from row-major to column-major for storage
      std::vector<std::vector<std::string>> columns_data(p, std::vector<std::string>(num_rows));
      
      for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < p; ++col) {
          columns_data[col][row] = matrix_data_rows[row][col];
        }
      }
      
      // Now push back the newly formatted columns
      for (size_t col = 0; col < p; ++col) {
        std::string col_name = base_name + "." + std::to_string(col + 1);
        
        if (containElementNamed(col_name)) 
          throw std::runtime_error("Column '" + col_name + "' already exists.");
        
        string_cols[col_name] = columns_data[col];
        names_.push_back(col_name);
      }
    }
  }
  
  
  void push_front(const std::vector<double>& col, const std::string& name) {
    if (containElementNamed(name))
      throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    numeric_cols[name] = col;
    names_.insert(names_.begin(), name);
  }
  
  void push_front(const std::vector<int>& col, const std::string& name) {
    if (containElementNamed(name))
      throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    int_cols[name] = col;
    names_.insert(names_.begin(), name);
  }
  
  void push_front(const std::vector<bool>& col, const std::string& name) {
    if (containElementNamed(name))
      throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    bool_cols[name] = col;
    names_.insert(names_.begin(), name);
  }
  
  void push_front(const std::vector<std::string>& col, const std::string& name) {
    if (containElementNamed(name))
      throw std::runtime_error("Column '" + name + "' already exists.");
    check_row_size(col.size(), name);
    string_cols[name] = col;
    names_.insert(names_.begin(), name);
  }
  
  /**
   * Pushes a double scalar to the front of the DataFrame by expanding it to a vector 
   * matching the DataFrame's current row count.
   */
  void push_front(double value, const std::string& name) {
    size_t current_rows = nrows();
    if (current_rows == 0 && !names_.empty()) {
      throw std::runtime_error("Cannot push scalar '" + name + "': DataFrame is currently empty (0 rows).");
    } else if (current_rows == 0) {
      // If DataFrame is empty, we can initialize it with one row
      current_rows = 1;
    }
    std::vector<double> expanded_col(current_rows, value);
    push_front(expanded_col, name); // Delegate to the vector overload
  }
  
  /**
   * Pushes an int scalar to the front of the DataFrame by expanding it to a vector 
   * matching the DataFrame's current row count.
   */
  void push_front(int value, const std::string& name) {
    size_t current_rows = nrows();
    if (current_rows == 0 && !names_.empty()) {
      throw std::runtime_error("Cannot push scalar '" + name + "': DataFrame is currently empty (0 rows).");
    } else if (current_rows == 0) {
      // If DataFrame is empty, we can initialize it with one row
      current_rows = 1;
    }
    std::vector<int> expanded_col(current_rows, value);
    push_front(expanded_col, name); // Delegate to the vector overload
  }
  
  /**
   * Pushes a bool scalar to the front of the DataFrame by expanding it to a vector 
   * matching the DataFrame's current row count.
   */
  void push_front(bool value, const std::string& name) {
    size_t current_rows = nrows();
    if (current_rows == 0 && !names_.empty()) {
      throw std::runtime_error("Cannot push scalar '" + name + "': DataFrame is currently empty (0 rows).");
    } else if (current_rows == 0) {
      // If DataFrame is empty, we can initialize it with one row
      current_rows = 1;
    }
    std::vector<bool> expanded_col(current_rows, value);
    push_front(expanded_col, name); // Delegate to the vector overload
  }
  
  /**
   * Pushes a string scalar to the front of the DataFrame by expanding it to a vector 
   * matching the DataFrame's current row count.
   */
  void push_front(const std::string& value, const std::string& name) {
    size_t current_rows = nrows();
    if (current_rows == 0 && !names_.empty()) {
      throw std::runtime_error("Cannot push scalar '" + name + "': DataFrame is currently empty (0 rows).");
    } else if (current_rows == 0) {
      // If DataFrame is empty, we can initialize it with one row
      current_rows = 1;
    }
    std::vector<std::string> expanded_col(current_rows, value);
    push_front(expanded_col, name); // Delegate to the vector overload
  }
  
  /**
   * Pushes a matrix (row-major std::vector<std::vector<double>>) to the front of the DataFrame.
   * Creates 1 column named 'base_name' if p==1, otherwise p columns named 'base_name'1 to 'base_name'p.
   */
  void push_front(const std::vector<std::vector<double>>& matrix_data_rows, const std::string& base_name) {
    size_t num_rows = matrix_data_rows.size();
    if (num_rows == 0) return;
    
    size_t p = matrix_data_rows[0].size();
    
    // Error checks and data transposition logic are the same as push_back
    for (size_t i = 1; i < num_rows; ++i) {
      if (matrix_data_rows[i].size() != p) {
        throw std::runtime_error("Input matrix has inconsistent column counts across rows.");
      }
    }
    check_row_size(num_rows, base_name);
    
    // We handle the insertion logic differently to use names_.insert()
    
    if (p == 1) {
      if (containElementNamed(base_name))
        throw std::runtime_error("Column '" + base_name + "' already exists.");
      
      std::vector<double> single_column(num_rows);
      for (size_t row = 0; row < num_rows; ++row) {
        single_column[row] = matrix_data_rows[row][0];
      }
      
      numeric_cols[base_name] = std::move(single_column);
      names_.insert(names_.begin(), base_name); // Insert name at the front
    }
    else {
      // p > 1: Transpose first
      std::vector<std::vector<double>> columns_data(p, std::vector<double>(num_rows));
      for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < p; ++col) {
          columns_data[col][row] = matrix_data_rows[row][col];
        }
      }
      
      // Insert names and data for multiple columns at the front
      std::vector<std::string> new_names;
      for (size_t col = 0; col < p; ++col) {
        std::string col_name = base_name + "." + std::to_string(col + 1);
        if (containElementNamed(col_name))
          throw std::runtime_error("Column '" + col_name + "' already exists.");
        
        numeric_cols[col_name] = columns_data[col];
        new_names.push_back(col_name);
      }
      // Insert all new names at the beginning of the names list
      names_.insert(names_.begin(), new_names.begin(), new_names.end());
    }
  }
  
  /**
   * Pushes a matrix (row-major std::vector<std::vector<int>>) to the front of the DataFrame.
   * Creates 1 column named 'base_name' if p==1, otherwise p columns named 'base_name'1 to 'base_name'p.
   */
  void push_front(const std::vector<std::vector<int>>& matrix_data_rows, const std::string& base_name) {
    size_t num_rows = matrix_data_rows.size();
    if (num_rows == 0) return;
    
    size_t p = matrix_data_rows[0].size();
    
    // Error checks and data transposition logic are the same as push_back
    for (size_t i = 1; i < num_rows; ++i) {
      if (matrix_data_rows[i].size() != p) {
        throw std::runtime_error("Input matrix has inconsistent column counts across rows.");
      }
    }
    check_row_size(num_rows, base_name);
    
    // We handle the insertion logic differently to use names_.insert()
    
    if (p == 1) {
      if (containElementNamed(base_name))
        throw std::runtime_error("Column '" + base_name + "' already exists.");
      
      std::vector<int> single_column(num_rows);
      for (size_t row = 0; row < num_rows; ++row) {
        single_column[row] = matrix_data_rows[row][0];
      }
      
      int_cols[base_name] = std::move(single_column);
      names_.insert(names_.begin(), base_name); // Insert name at the front
    }
    else {
      // p > 1: Transpose first
      std::vector<std::vector<int>> columns_data(p, std::vector<int>(num_rows));
      for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < p; ++col) {
          columns_data[col][row] = matrix_data_rows[row][col];
        }
      }
      
      // Insert names and data for multiple columns at the front
      std::vector<std::string> new_names;
      for (size_t col = 0; col < p; ++col) {
        std::string col_name = base_name + "." + std::to_string(col + 1);
        if (containElementNamed(col_name))
          throw std::runtime_error("Column '" + col_name + "' already exists.");
        
        int_cols[col_name] = columns_data[col];
        new_names.push_back(col_name);
      }
      // Insert all new names at the beginning of the names list
      names_.insert(names_.begin(), new_names.begin(), new_names.end());
    }
  }
  
  /**
   * Pushes a matrix (row-major std::vector<std::vector<bool>>) to the front of the DataFrame.
   * Creates 1 column named 'base_name' if p==1, otherwise p columns named 'base_name'1 to 'base_name'p.
   */
  void push_front(const std::vector<std::vector<bool>>& matrix_data_rows, const std::string& base_name) {
    size_t num_rows = matrix_data_rows.size();
    if (num_rows == 0) return;
    
    size_t p = matrix_data_rows[0].size();
    
    // Error checks and data transposition logic are the same as push_back
    for (size_t i = 1; i < num_rows; ++i) {
      if (matrix_data_rows[i].size() != p) {
        throw std::runtime_error("Input matrix has inconsistent column counts across rows.");
      }
    }
    check_row_size(num_rows, base_name);
    
    // We handle the insertion logic differently to use names_.insert()
    
    if (p == 1) {
      if (containElementNamed(base_name))
        throw std::runtime_error("Column '" + base_name + "' already exists.");
      
      std::vector<bool> single_column(num_rows);
      for (size_t row = 0; row < num_rows; ++row) {
        single_column[row] = matrix_data_rows[row][0];
      }
      
      bool_cols[base_name] = std::move(single_column);
      names_.insert(names_.begin(), base_name); // Insert name at the front
    }
    else {
      // p > 1: Transpose first
      std::vector<std::vector<bool>> columns_data(p, std::vector<bool>(num_rows));
      for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < p; ++col) {
          columns_data[col][row] = matrix_data_rows[row][col];
        }
      }
      
      // Insert names and data for multiple columns at the front
      std::vector<std::string> new_names;
      for (size_t col = 0; col < p; ++col) {
        std::string col_name = base_name + "." + std::to_string(col + 1);
        if (containElementNamed(col_name))
          throw std::runtime_error("Column '" + col_name + "' already exists.");
        
        bool_cols[col_name] = columns_data[col];
        new_names.push_back(col_name);
      }
      // Insert all new names at the beginning of the names list
      names_.insert(names_.begin(), new_names.begin(), new_names.end());
    }
  }
  
  /**
   * Pushes a matrix (row-major std::vector<std::vector<std::string>>) to the front of the DataFrame.
   * Creates 1 column named 'base_name' if p==1, otherwise p columns named 'base_name'1 to 'base_name'p.
   */
  void push_front(const std::vector<std::vector<std::string>>& matrix_data_rows, const std::string& base_name) {
    size_t num_rows = matrix_data_rows.size();
    if (num_rows == 0) return;
    
    size_t p = matrix_data_rows[0].size();
    
    // Error checks and data transposition logic are the same as push_back
    for (size_t i = 1; i < num_rows; ++i) {
      if (matrix_data_rows[i].size() != p) {
        throw std::runtime_error("Input matrix has inconsistent column counts across rows.");
      }
    }
    check_row_size(num_rows, base_name);
    
    // We handle the insertion logic differently to use names_.insert()
    
    if (p == 1) {
      if (containElementNamed(base_name))
        throw std::runtime_error("Column '" + base_name + "' already exists.");
      
      std::vector<std::string> single_column(num_rows);
      for (size_t row = 0; row < num_rows; ++row) {
        single_column[row] = matrix_data_rows[row][0];
      }
      
      string_cols[base_name] = std::move(single_column);
      names_.insert(names_.begin(), base_name); // Insert name at the front
    }
    else {
      // p > 1: Transpose first
      std::vector<std::vector<std::string>> columns_data(p, std::vector<std::string>(num_rows));
      for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < p; ++col) {
          columns_data[col][row] = matrix_data_rows[row][col];
        }
      }
      
      // Insert names and data for multiple columns at the front
      std::vector<std::string> new_names;
      for (size_t col = 0; col < p; ++col) {
        std::string col_name = base_name + "." + std::to_string(col + 1);
        if (containElementNamed(col_name))
          throw std::runtime_error("Column '" + col_name + "' already exists.");
        
        string_cols[col_name] = columns_data[col];
        new_names.push_back(col_name);
      }
      // Insert all new names at the beginning of the names list
      names_.insert(names_.begin(), new_names.begin(), new_names.end());
    }
  }
  
  
  void erase(const std::string& name) {
    if (numeric_cols.erase(name)) {
    } else if (int_cols.erase(name)) {
    } else if (bool_cols.erase(name)) {
    } else {
      string_cols.erase(name);
    }
    names_.erase(std::remove(names_.begin(), names_.end(), name), names_.end());
  }
  
  template <typename T>
  std::vector<T>& get(const std::string& name) {
    if constexpr (std::is_same_v<T, double>) {
      if (numeric_cols.count(name)) return numeric_cols.at(name);
    } else if constexpr (std::is_same_v<T, int>) {
      if (int_cols.count(name)) return int_cols.at(name);
    } else if constexpr (std::is_same_v<T, bool>) {
      if (bool_cols.count(name)) return bool_cols.at(name);
    } else if constexpr (std::is_same_v<T, std::string>) {
      if (string_cols.count(name)) return string_cols.at(name);
    }
    throw std::runtime_error("Column '" + name + "' not found or type mismatch.");
  }
  
  template <typename T>
  const std::vector<T>& get(const std::string& name) const {
    if constexpr (std::is_same_v<T, double>) {
      if (numeric_cols.count(name)) return numeric_cols.at(name);
    } else if constexpr (std::is_same_v<T, int>) {
      if (int_cols.count(name)) return int_cols.at(name);
    } else if constexpr (std::is_same_v<T, bool>) {
      if (bool_cols.count(name)) return bool_cols.at(name);
    } else if constexpr (std::is_same_v<T, std::string>) {
      if (string_cols.count(name)) return string_cols.at(name);
    }
    throw std::runtime_error("Column '" + name + "' not found or type mismatch.");
  }
};



// Define ListCpp to hold heterogeneous data using std::unordered_map and std::variant
struct ListCpp {
  std::vector<std::string> names_;
  std::unordered_map<std::string, std::variant<
    bool, int, double, std::string,
    std::vector<bool>,
    std::vector<int>,
    std::vector<double>,
    std::vector<std::string>,
    std::vector<std::vector<bool>>,
    std::vector<std::vector<int>>,
    std::vector<std::vector<double>>,
    std::vector<std::vector<std::string>>,
    DataFrameCpp,
    ListCpp,
    std::vector<DataFrameCpp>,
    std::vector<ListCpp>
    >> data;
    
    // Return the number of elements in the list
    size_t size() const {
      return data.size();
    }
    
    // Check if an element exists by name
    bool containsElementNamed(const std::string& name) const {
      return data.count(name) > 0;
    }
    
    // Push an item to the end, preserving order in names_.
    template<typename T>
    void push_back(T value, const std::string& name) {
      if (containsElementNamed(name))
        throw std::runtime_error("Column '" + name + "' already exists.");
      data[name] = value;
      names_.push_back(name);
    }
    
    // Insert an item at the beginning, preserving order in names_.
    template<typename T>
    void push_front(T value, const std::string& name) {
      if (containsElementNamed(name))
        throw std::runtime_error("Column '" + name + "' already exists.");
      data[name] = value;
      names_.insert(names_.begin(), name);
    }
    
    // Return the names (keys) of the list elements in insertion order
    std::vector<std::string> names() const {
      return names_;
    }
    
    // Get an element by name and type
    template <typename T>
    T& get(const std::string& name) {
      if (!containsElementNamed(name))
        throw std::runtime_error("Element with name '" + name + "' not found.");
      return std::get<T>(data.at(name));
    }
    
    // Get an element by name and type (const version)
    template <typename T>
    const T& get(const std::string& name) const {
      if (!containsElementNamed(name))
        throw std::runtime_error("Element with name '" + name + "' not found.");
      return std::get<T>(data.at(name));
    }
    
    // Erase an element by name
    void erase(const std::string& name) {
      data.erase(name);
      names_.erase(std::remove(names_.begin(), names_.end(), name), names_.end());
    }
    
};


// Convert Rcpp::DataFrame to DataFrameCpp
inline DataFrameCpp convertRDataFrameToCpp(const Rcpp::DataFrame& r_df) {
  DataFrameCpp df_cpp;
  Rcpp::CharacterVector names = r_df.names();
  size_t num_cols = r_df.size();
  
  if (num_cols == 0) {
    return df_cpp;
  }
  
  // Iterate through all columns in the Rcpp DataFrame
  for (size_t i = 0; i < num_cols; ++i) {
    std::string col_name = Rcpp::as<std::string>(names[i]);
    Rcpp::RObject current_col = r_df[i];
    
    // Check the type of the R column and use the appropriate push_back overload
    if (Rcpp::is<Rcpp::NumericVector>(current_col)) {
      // Rcpp::NumericVector wraps std::vector<double> naturally
      Rcpp::NumericVector vec = Rcpp::as<Rcpp::NumericVector>(current_col);
      df_cpp.push_back(Rcpp::as<std::vector<double>>(vec), col_name);
      
    } else if (Rcpp::is<Rcpp::IntegerVector>(current_col)) {
      // Rcpp::IntegerVector wraps std::vector<int> naturally
      Rcpp::IntegerVector vec = Rcpp::as<Rcpp::IntegerVector>(current_col);
      df_cpp.push_back(Rcpp::as<std::vector<int>>(vec), col_name);
      
    } else if (Rcpp::is<Rcpp::LogicalVector>(current_col)) {
      // Rcpp::LogicalVector wraps std::vector<bool> naturally
      Rcpp::LogicalVector vec = Rcpp::as<Rcpp::LogicalVector>(current_col);
      df_cpp.push_back(Rcpp::as<std::vector<bool>>(vec), col_name);
      
    } else if (Rcpp::is<Rcpp::CharacterVector>(current_col)) {
      // Rcpp::CharacterVector wraps std::vector<std::string> naturally
      Rcpp::CharacterVector vec = Rcpp::as<Rcpp::CharacterVector>(current_col);
      df_cpp.push_back(Rcpp::as<std::vector<std::string>>(vec), col_name);
      
    } else {
      // Handle unsupported types if necessary (e.g., factors, lists within DF)
      Rcpp::warning("Column type not supported in DataFrameCpp conversion: " + col_name);
    }
  }
  
  return df_cpp;
}



// Forward declarations for the conversion functions (needed for recursion)
inline Rcpp::List convertListCppToR(const ListCpp& list_cpp);

// Convert DataFrameCpp to Rcpp::DataFrame
inline Rcpp::DataFrame convertDataFrameCppToR(const DataFrameCpp& df_cpp) {
  // Rcpp::List will build the underlying structure of the R data frame
  Rcpp::List columns;
  Rcpp::CharacterVector col_names;
  
  // Iterate through all column names in the order they were inserted
  for (const std::string& name : df_cpp.names_) {
    col_names.push_back(name);
    
    // Check which map the column belongs to and add it to the Rcpp::List
    if (df_cpp.numeric_cols.count(name)) {
      // Rcpp::NumericVector implicitly converts std::vector<double>
      columns.push_back(Rcpp::wrap(df_cpp.numeric_cols.at(name)));
    } else if (df_cpp.int_cols.count(name)) {
      // Rcpp::IntegerVector implicitly converts std::vector<int>
      columns.push_back(Rcpp::wrap(df_cpp.int_cols.at(name)));
    } else if (df_cpp.bool_cols.count(name)) {
      // Rcpp::LogicalVector implicitly converts std::vector<bool>
      columns.push_back(Rcpp::wrap(df_cpp.bool_cols.at(name)));
    } else if (df_cpp.string_cols.count(name)) {
      // Rcpp::CharacterVector implicitly converts std::vector<std::string>
      columns.push_back(Rcpp::wrap(df_cpp.string_cols.at(name)));
    }
  }
  
  // Assign names to the list of columns
  columns.names() = col_names;
  
  // Create the Rcpp::DataFrame from the list of columns and set row names/classes
  // Assuming all columns have the same length (which your check_row_size enforces)
  return Rcpp::as<Rcpp::DataFrame>(columns);
}



// Updated RcppVisitor struct:
struct RcppVisitor {
  // Make every function explicitly return Rcpp::RObject
  
  // Overload for simple types and std::vectors (handled by Rcpp::wrap)
  template<typename T>
  Rcpp::RObject operator()(const T& value) const {
    return Rcpp::wrap(value);
  }

  // Overload for DataFrameCpp (explicitly return Rcpp::RObject)
  Rcpp::RObject operator()(const DataFrameCpp& df) const {
    return convertDataFrameCppToR(df); 
  }
  
  // Overload for std::vector<DataFrameCpp> 
  Rcpp::RObject operator()(const std::vector<DataFrameCpp>& dfs) const {
    Rcpp::List result;
    for (const auto& df : dfs) {
      result.push_back(convertDataFrameCppToR(df));
    }
    return result; 
  }
  
  // Overload for recursive ListCpp
  Rcpp::RObject operator()(const ListCpp& list) const {
    return convertListCppToR(list); 
  }
  
  // Overload for std::vector<ListCpp>
  Rcpp::RObject operator()(const std::vector<ListCpp>& lists) const {
    Rcpp::List result;
    for (const auto& list : lists) {
      result.push_back(convertListCppToR(list));
    }
    return result; 
  }
};



// Convert ListCpp to Rcpp::List
inline Rcpp::List convertListCppToR(const ListCpp& list_cpp) {
  Rcpp::List r_list;
  
  // Ensure we maintain the insertion order using the names_ vector
  for (const std::string& name : list_cpp.names_) {
    // Get the variant object from the unordered map
    const auto& variant_value = list_cpp.data.at(name);
    
    // Use std::visit to apply the correct conversion operator from the visitor struct
    Rcpp::RObject r_object = std::visit(RcppVisitor{}, variant_value);
    
    r_list.push_back(r_object, name); // Add the converted object to the R list with the correct name
  }
  
  return r_list;
}



namespace Rcpp {
template <> inline SEXP wrap(const DataFrameCpp& df_cpp) {
  // Now the compiler knows what both types and the function are.
  return Rcpp::wrap(convertDataFrameCppToR(df_cpp));
}

template <> inline SEXP wrap(const ListCpp& list_cpp) {
  return Rcpp::wrap(convertListCppToR(list_cpp));
}
}


#endif // __DATAFRAME_LIST__