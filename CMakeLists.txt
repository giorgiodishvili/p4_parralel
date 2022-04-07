cmake_minimum_required(VERSION 3.21)
project(p4_parallel C)

#set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "") # works
set(CMAKE_C_STANDARD 99)

add_executable(p4_parallel main.c main2-a.c main2-b.c)

#find_package(MPI REQUIRED PATHS /opt/homebrew/Cellar/open-mpi/4.1.2)
include_directories(/opt/homebrew/Cellar/open-mpi/4.1.2/include)