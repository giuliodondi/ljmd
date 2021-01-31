#include "gtest/gtest.h"


#include "physics.h"





class Ekin_T_Test : public ::testing::Test {

protected:
        mdsys_t *sys;

        void SetUp() {
                sys = new mdsys_t;
                sys->natoms = 2;
                sys->mass = 1.0;
                sys->vx = new double[2];
				sys->vy = new double[2];
				sys->vz = new double[2];

                sys->vx[0] = 0.0;
                sys->vx[1] = 0.0;
				sys->vy[0] = 1.0;
                sys->vy[1] = -1.0;
				sys->vz[0] = 0.0;
                sys->vz[1] = 0.0;

        }

        void TearDown() {
                delete[] sys->vx;
				delete[] sys->vy;
				delete[] sys->vz;

                delete sys;
        }
};

TEST(Ekin_T_TestEmpty, empty) {
        mdsys_t *sys = new mdsys_t;
        sys->natoms = 0;

        ekin(sys);
	
        ASSERT_DOUBLE_EQ(sys->ekin, 0.0);
        ASSERT_DOUBLE_EQ(sys->temp, 0.0);

        delete sys;
}



TEST_F(Ekin_T_Test, test1) {
        ASSERT_NE(sys, nullptr);
        ASSERT_DOUBLE_EQ(sys->vy[0], 1.0);
        ASSERT_DOUBLE_EQ(sys->vy[1], -1.0);
	
        ekin(sys);
		double exp_ekin = mvsq2e;
		double exp_temp = 2*exp_ekin / (3* kboltz );
	
        ASSERT_DOUBLE_EQ(sys->ekin, exp_ekin);
        ASSERT_DOUBLE_EQ(sys->temp, exp_temp);
}
