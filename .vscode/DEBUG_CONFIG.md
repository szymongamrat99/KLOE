# Debug Configuration Guide - KLSPM00 Project

## Overview

Configured debug launch configurations for C/C++ development with CERN ROOT and Boost libraries.

## Available Debug Configurations

### 1. **KLSPM00 Debug (GDB)** - Recommended for standard debugging
- **Purpose**: Basic debugging with GDB
- **Features**:
  - Exception breakpoint on throw
  - Pretty-printing enabled
  - Address Sanitizer enabled
  - Auto-builds in Debug mode before launch

### 2. **KLSPM00 Debug (GDB - ROOT pretty-print)**
- **Purpose**: Enhanced debugging with ROOT object pretty-printers
- **Features**:
  - Loads ROOT's gdb-print.py for better ROOT object inspection
  - Exception breakpoint on throw
  - Useful for inspecting TH1F, TTree, and other ROOT objects
  - Auto-builds in Debug mode before launch

### 3. **KLSPM00 Debug with AddressSanitizer**
- **Purpose**: Memory leak and undefined behavior detection
- **Features**:
  - AddressSanitizer enabled for leak detection
  - Breaks on sanitizer errors
  - Shows detailed memory error information
  - Great for catching memory bugs early

### 4. **KLSPM00 Release Debug**
- **Purpose**: Debug optimized Release build
- **Features**:
  - Release build with optimization (-O3)
  - Debugging symbols included
  - Use when you need performance with some debugging capability

## How to Use

### Starting Debug Session
1. Press `F5` or go to **Run and Debug** in the left sidebar
2. Select desired configuration from the dropdown
3. The build task will run automatically before launching
4. Debugger will start and stop at breakpoints

### Setting Breakpoints
- Click on the line number margin to set/remove breakpoints
- Use conditional breakpoints: Right-click → Add Conditional Breakpoint
- Use function name breakpoints in Debug Console: `break function_name`

### Debug Console Commands
```
# Common GDB commands in Debug Console
break filename.cpp:123          # Set breakpoint at line
watch variable_name             # Watch variable changes
continue/cont                   # Continue execution
step                           # Step into function
next                           # Step over
print variable_name            # Print variable value
print *this                    # Print object members
backtrace                      # Show call stack
```

### Inspecting ROOT Objects
When using "ROOT pretty-print" configuration:
```
print my_histogram->GetEntries()
print my_tree->GetEntries()
print my_roofit_var
```

## Build Tasks

Available tasks in Command Palette (`Ctrl+Shift+P`):
- **build-debug**: Build with Debug flags and sanitizers
- **build-release**: Build with Release optimization
- **rebuild-debug**: Clean and rebuild Debug build
- **clean**: Clean build artifacts

## Environment Configuration

### LD_LIBRARY_PATH
Points to:
- CERN ROOT: `/data/4/users/gamrat/root_dirs/install_root/lib`
- System libraries: `/usr/lib/x86_64-linux-gnu`
- Custom: `$LD_LIBRARY_PATH`

### AddressSanitizer Options
```
ASAN_OPTIONS=abort_on_error=1:detect_leaks=1:verbosity=2
```

## Troubleshooting

### Library Loading Issues
If you get "libCore.so not found":
1. Ensure LD_LIBRARY_PATH is set correctly
2. Verify ROOT installation at: `/data/4/users/gamrat/root_dirs/install_root/lib`
3. Check: `ldd ./build/bin/KLSPM00`

### GDB Not Found
Install GDB:
```bash
sudo apt-get install gdb
```

### Poor Symbol Resolution
1. Ensure build type is Debug: `-g3 -O0`
2. Rebuild with `rebuild-debug` task
3. Check: `nm ./build/bin/KLSPM00 | grep function_name`

### AddressSanitizer Linking Error
Ensure CMakeLists.txt has:
```cmake
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -fsanitize=address -fsanitize=undefined")
```

## Tips & Best Practices

1. **Use Debug mode by default** for development
2. **Enable pretty-printers** for ROOT projects
3. **Set conditional breakpoints** to avoid stopping too often
4. **Use Watch tab** to monitor key variables
5. **Check Memory leaks** periodically with AddressSanitizer config
6. **Debug with Release build** only if needed for performance analysis

## File Locations

- **Launch configs**: `.vscode/launch.json`
- **Tasks configs**: `.vscode/tasks.json`
- **IntelliSense config**: `.vscode/c_cpp_properties.json`
- **Executable**: `build/bin/KLSPM00`
- **ROOT installation**: `/data/4/users/gamrat/root_dirs/install_root/`
- **Boost libraries**: `/usr/lib/x86_64-linux-gnu/`

## Additional Resources

- [VS Code Debugging](https://code.visualstudio.com/docs/editor/debugging)
- [GDB Documentation](https://sourceware.org/gdb/documentation/)
- [CERN ROOT Documentation](https://root.cern/manual/)
- [AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer)
