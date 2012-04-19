################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../code/function.o \
../code/main.o \
../code/solver.o \
../code/thomas.o 

CPP_SRCS += \
../code/function.cpp \
../code/main.cpp \
../code/solver.cpp \
../code/thomas.cpp 

CC_SRCS += \
../code/solve.cc 

OBJS += \
./code/function.o \
./code/main.o \
./code/solve.o \
./code/solver.o \
./code/thomas.o 

CC_DEPS += \
./code/solve.d 

CPP_DEPS += \
./code/function.d \
./code/main.d \
./code/solver.d \
./code/thomas.d 


# Each subdirectory must supply rules for building sources it contributes
code/%.o: ../code/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

code/%.o: ../code/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


