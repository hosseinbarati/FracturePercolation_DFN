//#include "GPU_SOLVER.h"
//
//
//
//GPU_SOLVER::GPU_SOLVER()
//{
//	cl::Platform::get(&platforms);
//	std::cout << "There are " << platforms.size() << " platforms on your machine." << std::endl;
//
//	for (size_t i = 0; i < platforms.size(); i++) {
//		std::cout << "Platform " << i << endl;
//		std::cout << "CL_PLATFORM_VENDOR: " << platforms[i].getInfo<CL_PLATFORM_VENDOR>() << std::endl;
//		std::cout << "CL_PLATFORM_NAME: " << platforms[i].getInfo<CL_PLATFORM_NAME>() << std::endl;
//		std::cout << "CL_PLATFORM_VERSION: " << platforms[i].getInfo<CL_PLATFORM_VERSION>() << std::endl;
//		std::cout << "CL_PLATFORM_PROFILE: " << platforms[i].getInfo<CL_PLATFORM_PROFILE>() << std::endl;
//		std::cout << "CL_PLATFORM_EXTENSIONS: " << platforms[i].getInfo<CL_PLATFORM_EXTENSIONS>() << std::endl;
//		std::cout << std::endl;
//		platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);
//
//		for (size_t j = 0; j < devices.size(); j++) {
//			std::cout << "Platform " << i << "		Devices " << j << std::endl;
//			std::cout << "CL_DEVICE_NAME: " << devices[j].getInfo<CL_DEVICE_NAME>() << std::endl;
//			std::cout << "CL_DEVICE_VENDOR: " << devices[j].getInfo<CL_DEVICE_VENDOR>() << std::endl;
//			std::cout << "CL_DEVICE_MAX_COMPUTE_UNITS: " << devices[j].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
//			std::cout << "CL_DEVICE_MAX_CLOCK_FREQUENCY: " << devices[j].getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
//			std::cout << "CL_DEVICE_LOCAL_MEM_SIZE: " << devices[j].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
//			std::cout << "CL_DEVICE_GLOBAL_MEM_SIZE: " << devices[j].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << std::endl;
//		}
//	}
//
//	
//
//
//	
//}
//
//
//void GPU_SOLVER::creatContextFromType(cl_device_type devType, int i) {
//
//	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[i])(), 0 };
//	context = cl::Context(devType, cps);
//}
//
//void GPU_SOLVER::creatContextFromIndex(int pidx, int didx) {
//	std::vector < cl::Device > device;
//	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[pidx])(), 0 };
//	device.push_back(devices[didx]);
//	context = cl::Context(device, cps, NULL, NULL);
//}
//
//void GPU_SOLVER::listDevicesInContext() {
//	std::vector<cl::Device> tempDevices = context.getInfo<CL_CONTEXT_DEVICES>();
//	for (size_t i = 0; i < devices.size(); i++) {
//		std::cout << "Device " << i << std::endl;
//		std::cout << "CL_DEVICE_NAME: " << devices[i].getInfo <CL_DEVICE_NAME> () << std::endl;
//		std::cout << "CL_DEVICE_VENDOR: " << devices[i].getInfo < CL_DEVICE_VENDOR > () << std::endl;
//	}
//}
//
//
//GPU_SOLVER::~GPU_SOLVER()
//{
//}
