################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/benchmarks/rings.cpp \
../src/benchmarks/rotation.cpp \
../src/benchmarks/tensile.cpp 

OBJS += \
./src/benchmarks/rings.o \
./src/benchmarks/rotation.o \
./src/benchmarks/tensile.o 

CPP_DEPS += \
./src/benchmarks/rings.d \
./src/benchmarks/rotation.d \
./src/benchmarks/tensile.d 


# Each subdirectory must supply rules for building sources it contributes
src/benchmarks/%.o: ../src/benchmarks/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++14 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


