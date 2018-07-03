################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/element.cpp \
../src/gauss_integration.cpp \
../src/kd_tree.cpp \
../src/kernel.cpp \
../src/mfree_iwf.cpp \
../src/particle.cpp \
../src/simulation_data.cpp \
../src/simulation_time.cpp 

OBJS += \
./src/element.o \
./src/gauss_integration.o \
./src/kd_tree.o \
./src/kernel.o \
./src/mfree_iwf.o \
./src/particle.o \
./src/simulation_data.o \
./src/simulation_time.o 

CPP_DEPS += \
./src/element.d \
./src/gauss_integration.d \
./src/kd_tree.d \
./src/kernel.d \
./src/mfree_iwf.d \
./src/particle.d \
./src/simulation_data.d \
./src/simulation_time.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++14 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


