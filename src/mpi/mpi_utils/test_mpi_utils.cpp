#include "gtest/gtest.h"
#include <math.h>

#include "mpi_headers/mpi_utils.h"

TEST(right_triag_area_test, max_position) {

        int l = 9;
        int area = right_triangle_area(l);
        EXPECT_EQ(area, 45);
}

class max_min_test : public ::testing::TestWithParam<int> {

      protected:
        int max_idx = GetParam();
};

TEST_P(max_min_test, max_position) {

        int array[4] = {0, 0, 0, 0};

        array[max_idx] = 2;

        int max, min;

        max_min_arr(array, 4, &max, &min);

        EXPECT_EQ(min, 0);
        EXPECT_EQ(max, 2);
}

INSTANTIATE_TEST_SUITE_P(max_min_test_parametric, max_min_test,
                         ::testing::Values(0, 1, 2, 3));

class maxmin_idx_test : public ::testing::TestWithParam<std::tuple<int, int>> {

      protected:
        int val = std::get<0>(GetParam());
        int val_idx = std::get<1>(GetParam());
};

TEST_P(maxmin_idx_test, both_cycle) {

        int array[4] = {0, 0, 0, 0};

        array[val_idx] = val;

        int swch = 0;

        if (val > 0) {
                swch = 1;
        }

        int out_idx = max_min_index(array, 4, swch);

        EXPECT_EQ(out_idx, val_idx);
}

INSTANTIATE_TEST_SUITE_P(maxmin_idx_test_parametric, maxmin_idx_test,
                         ::testing::Combine(::testing::Values(2, -2),
                                            ::testing::Values(0, 1, 2, 3)));

TEST(split_triangle, single_test) {

        int nprocs = 3;
        int size = 8;

        int segment[3 * nprocs];
        int segment_areas[3 * nprocs];

        split_triangle_equal_areas(size, nprocs, segment, segment_areas);

        EXPECT_EQ(segment[0], 4);
        EXPECT_EQ(segment[1], 6);
        EXPECT_EQ(segment[2], 8);
        EXPECT_EQ(segment[3], 5);
        EXPECT_EQ(segment[4], 7);
        EXPECT_EQ(segment[5], 8);
        EXPECT_EQ(segment[6], 5);
        EXPECT_EQ(segment[7], 7);
        EXPECT_EQ(segment[8], 8);

        EXPECT_EQ(segment_areas[0], 10);
        EXPECT_EQ(segment_areas[1], 11);
        EXPECT_EQ(segment_areas[2], 15);
        EXPECT_EQ(segment_areas[3], 15);
        EXPECT_EQ(segment_areas[4], 13);
        EXPECT_EQ(segment_areas[5], 8);
        EXPECT_EQ(segment_areas[6], 15);
        EXPECT_EQ(segment_areas[7], 13);
        EXPECT_EQ(segment_areas[8], 8);
}

class segments_test : public ::testing::TestWithParam<int> {

      protected:
        arr_seg_t *proc_seg;
        int nprocs;
        int size;
        int proc_id = GetParam();

        void SetUp() {
                nprocs = 3;
                size = 8;
                proc_seg = new arr_seg_t;
                proc_seg->splitting = new int[nprocs]();
        }

        void TearDown() {

                delete proc_seg->splitting;
                delete proc_seg;
        }
};

TEST_P(segments_test, init) {

        ASSERT_EQ(nprocs, 3);
        ASSERT_EQ(size, 8);
        ASSERT_TRUE(proc_id >= 0);
        ASSERT_TRUE(proc_id < nprocs);
        ASSERT_TRUE(proc_seg);
        ASSERT_TRUE(proc_seg);
        ASSERT_TRUE(proc_seg->splitting);

        init_segments(nprocs, proc_id, proc_seg, size);

        EXPECT_EQ(proc_seg->splitting[0], 2);
        EXPECT_EQ(proc_seg->splitting[1], 2);
        EXPECT_EQ(proc_seg->splitting[2], 4);

        EXPECT_EQ(proc_seg->size, proc_seg->splitting[proc_id]);

        int exp_idx = 0;
        for (int p = 0; p < proc_id; ++p) {
                exp_idx += proc_seg->splitting[p];
        }
        EXPECT_EQ(proc_seg->idx, exp_idx);
}

TEST_P(segments_test, comm_arrays) {

        ASSERT_EQ(nprocs, 3);
        ASSERT_EQ(size, 8);
        ASSERT_TRUE(proc_id >= 0);
        ASSERT_TRUE(proc_id < nprocs);
        ASSERT_TRUE(proc_seg);
        ASSERT_TRUE(proc_seg);
        ASSERT_TRUE(proc_seg->splitting);

        init_segments(nprocs, proc_id, proc_seg, size);

        int count[nprocs];
        int offsets[nprocs];

        mpi_collective_comm_arrays(nprocs, proc_seg->splitting, count, offsets);

        EXPECT_EQ(count[0], 2);
        EXPECT_EQ(count[1], 2);
        EXPECT_EQ(count[2], 4);

        EXPECT_EQ(offsets[0], 0);
        EXPECT_EQ(offsets[1], 2);
        EXPECT_EQ(offsets[2], 4);
}

INSTANTIATE_TEST_SUITE_P(segments_test_parametric, segments_test,
                         ::testing::Values(0, 1, 2));