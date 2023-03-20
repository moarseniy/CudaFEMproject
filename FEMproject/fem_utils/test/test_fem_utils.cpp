#include <fem_utils/fem_utils.h>
#include <cuda_fem_utils/cuda_fem_utils.h>

#include <gtest/gtest.h>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
#include <unordered_set>
#include <map>

TEST(UtilsFuncs, test) {
  EXPECT_EQ(true, true);
}

int main(int argc, char *argv[]) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
