#include <TestMain.hpp>

Hello::Hello() {
	std::cout << "hi!" << std::endl;
}

Hello::Hello(std::string s) {
	std::cout << "hi " << s << std::endl;
}

int Hello::get() {
	return 42;
}