#include "gtest/gtest.h"


extern "C" {
#include "integration.c"
}


class IntegrationTest : public ::testing::TestWithParam<std::tuple<double,double>> {

protected:
        mdsys_t *sys;
		//double force_params = GetParam();

        void SetUp() {
                sys = new mdsys_t;
                sys->natoms = 2;
                sys->mass = 1.0;
			
                sys->rx = new double[2]();
                sys->ry = new double[2]();
                sys->rz = new double[2]();
			
				sys->vx = new double[2]();
                sys->vy = new double[2]();
                sys->vz = new double[2]();

                sys->fx = new double[2]();
                sys->fy = new double[2]();
                sys->fz = new double[2]();
			
        }

        void TearDown() {
                delete[] sys->rx;
                delete[] sys->ry;
                delete[] sys->rz;
			
				delete[] sys->vx;
                delete[] sys->vy;
                delete[] sys->vz;

                delete[] sys->fx;
                delete[] sys->fy;
                delete[] sys->fz;

                delete sys;
        }
};

TEST_P(IntegrationTest, testVerlet1) {
        ASSERT_NE(sys, nullptr);
	
		sys->dt = std::get<0>(GetParam());
	
        sys->ry[0] = 1.0;
		sys->ry[1] = -1.0;
	
		sys->fx[0] = std::get<1>(GetParam());
		sys->fx[1] = std::get<1>(GetParam());
	
		verlet_1(sys);
	
		double exp_v_increm = 0.5*std::get<0>(GetParam())*std::get<1>(GetParam()) / mvsq2e;
		double exp_r_increm = exp_v_increm*std::get<0>(GetParam());
		
		ASSERT_DOUBLE_EQ(sys->rx[0],exp_r_increm);
		ASSERT_DOUBLE_EQ(sys->ry[0],1.0);
		ASSERT_DOUBLE_EQ(sys->rz[0],0.0);
		ASSERT_DOUBLE_EQ(sys->vx[0],exp_v_increm);
		ASSERT_DOUBLE_EQ(sys->vy[0],0.0);
		ASSERT_DOUBLE_EQ(sys->vz[0],0.0);

		ASSERT_DOUBLE_EQ(sys->rx[1],exp_r_increm);
		ASSERT_DOUBLE_EQ(sys->ry[1],-1.0);
		ASSERT_DOUBLE_EQ(sys->rz[1],0.0);
		ASSERT_DOUBLE_EQ(sys->vx[1],exp_v_increm);
		ASSERT_DOUBLE_EQ(sys->vy[1],0.0);
		ASSERT_DOUBLE_EQ(sys->vz[1],0.0);

}

TEST_P(IntegrationTest, testVerlet2) {
        ASSERT_NE(sys, nullptr);
	
		sys->dt = std::get<0>(GetParam());
	
        sys->ry[0] = 1.0;
		sys->ry[1] = -1.0;
	
		sys->fx[0] = std::get<1>(GetParam());
		sys->fx[1] = std::get<1>(GetParam());
	
		verlet_2(sys);
	
		double exp_v_increm = 0.5*std::get<0>(GetParam())*std::get<1>(GetParam()) / mvsq2e;
		
		ASSERT_DOUBLE_EQ(sys->rx[0],0.0);
		ASSERT_DOUBLE_EQ(sys->ry[0],1.0);
		ASSERT_DOUBLE_EQ(sys->rz[0],0.0);
		ASSERT_DOUBLE_EQ(sys->vx[0],exp_v_increm);
		ASSERT_DOUBLE_EQ(sys->vy[0],0.0);
		ASSERT_DOUBLE_EQ(sys->vz[0],0.0);

		ASSERT_DOUBLE_EQ(sys->rx[1],0.0);
		ASSERT_DOUBLE_EQ(sys->ry[1],-1.0);
		ASSERT_DOUBLE_EQ(sys->rz[1],0.0);
		ASSERT_DOUBLE_EQ(sys->vx[1],exp_v_increm);
		ASSERT_DOUBLE_EQ(sys->vy[1],0.0);
		ASSERT_DOUBLE_EQ(sys->vz[1],0.0);

}

INSTANTIATE_TEST_SUITE_P(IntegrationTest_parametric,
			 IntegrationTest,
			 ::testing::Combine(
				::testing::Values(0.0, 0.5, 1.0),
			 	::testing::Values(0.0, 1.0, 2.0))
						);