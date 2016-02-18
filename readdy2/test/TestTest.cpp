#include <TestMain.hpp>
#include "gtest/gtest.h"

namespace{
	class TestTest : public ::testing::Test {
	protected:
		Hello *h;

	  virtual void SetUp() {
		  h = new Hello();
		  std::cout << "set up test!" << std::endl;
	  }

	  virtual void TearDown() {
		  delete h;
		  std::cout << "tear down test!" << std::endl;
	  }

	};
	
	TEST_F(TestTest, ReturnsTheRightThing) {
		std::cout << "testing 42" << std::endl;
		EXPECT_EQ(42, h->get());
	}
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}