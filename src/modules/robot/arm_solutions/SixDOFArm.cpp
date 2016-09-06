#include "SixDOFArmSolution.h"
#include <fastmath.h>
#include "checksumm.h"
#include "ActuatorCoordinates.h"
#include "ConfigValue.h"
#include "libs/Kernel.h"
#include "StreamOutputPool.h"

#include "libs/nuts_bolts.h"

#include "libs/Config.h"

#define l1_length_checksum          CHECKSUM("l1_length")
#define l2_length_checksum          CHECKSUM("l2_length")
#define l3_length_checksum          CHECKSUM("l3_length")
#define l4_length_checksum          CHECKSUM("l4_length")
#define homing_checksum        CHECKSUM("homing")

#define SQ(x) powf(x, 2)
#define ROUND(x, y) (roundf(x * 1e ## y) / 1e ## y)

SixDOFArmSolution::SixDOFArmSolution(Config* config)
{
    // arm1_length is the length of the inner main arm from hinge to hinge
    l1_length         = config->value(l1_length_checksum)->by_default(202.0f)->as_number();
    // arm2_length is the length of the inner main arm from hinge to hinge
    l2_length         = config->value(l2_length_checksum)->by_default(160.0f)->as_number();
    // arm1_length is the length of the inner main arm from hinge to hinge
    l3_length         = config->value(l3_length_checksum)->by_default(195.0f)->as_number();
    // arm1_length is the length of the inner main arm from hinge to hinge
    l4_length         = config->value(l4_length_checksum)->by_default(67.15f)->as_number();

    init();
}

void SixDOFArmSolution::init()
{

}

float SixDOFArmSolution::to_degrees(float radians) const
{
    return radians * (180.0F / 3.14159265359f);
}

void SixDOFArmSolution::cartesian_to_actuator(const float cartesian_mm[], ActuatorCoordinates &actuator_mm ) const
{

    float tcp_pos[3],
          aux_pos[3],
          tcp_ori[3],
          ARM_Q1,
          ARM_Q2,
          ARM_Q3,
          ARM_Q4,
          ARM_Q5,
          ARM_Q6,
          C1,
          C2,
          C3,
          S1,
          S2,
          S3;

    tcp_pos[X_AXIS] = cartesian_mm[X_AXIS];
    tcp_pos[Y_AXIS] = cartesian_mm[Y_AXIS];
    tcp_pos[Z_AXIS] = cartesian_mm[Z_AXIS];

    /*aux_pos[X_AXIS] = tcp_pos[X_AXIS]-(this->l4_length*tcp_ori[X_AXIS]);
    aux_pos[Y_AXIS] = tcp_pos[Y_AXIS]-(this->l4_length*tcp_ori[Y_AXIS]);
    aux_pos[Z_AXIS] = tcp_pos[Z_AXIS]-(this->l4_length*tcp_ori[Z_AXIS]);*/

    aux_pos[X_AXIS] = tcp_pos[X_AXIS];
    aux_pos[Y_AXIS] = tcp_pos[Y_AXIS];  //TEMPORAL
    aux_pos[Z_AXIS] = tcp_pos[Z_AXIS];

    /*
    if (this->arm1_length == this->arm2_length)
        SCARA_C2 = (SQ(SCARA_pos[X_AXIS]) + SQ(SCARA_pos[Y_AXIS]) - 2.0f * SQ(this->arm1_length)) / (2.0f * SQ(this->arm1_length));
    else
        SCARA_C2 = (SQ(SCARA_pos[X_AXIS]) + SQ(SCARA_pos[Y_AXIS]) - SQ(this->arm1_length) - SQ(this->arm2_length)) / (2.0f * SQ(this->arm1_length));

    // SCARA position is undefined if abs(SCARA_C2) >=1
    // In reality abs(SCARA_C2) >0.95 can be problematic.

    if (SCARA_C2 > this->morgan_undefined_max)
        SCARA_C2 = this->morgan_undefined_max;
    else if (SCARA_C2 < -this->morgan_undefined_min)
        SCARA_C2 = -this->morgan_undefined_min;


    SCARA_S2 = sqrtf(1.0f - SQ(SCARA_C2));

    SCARA_K1 = this->arm1_length + this->arm2_length * SCARA_C2;
    SCARA_K2 = this->arm2_length * SCARA_S2;

    SCARA_theta = (atan2f(SCARA_pos[X_AXIS], SCARA_pos[Y_AXIS]) - atan2f(SCARA_K1, SCARA_K2)) * -1.0f; // Morgan Thomas turns Theta in oposite direction
    SCARA_psi   = atan2f(SCARA_S2, SCARA_C2);*/

    ARM_Q1=atan(aux_pos[X_AXIS]/aux_pos[Y_AXIS]);
    C1=cosf(ARM_Q1);
    S1=sinf(ARM_Q1);
    C3=acos((SQ(aux_pos[X_AXIS])+SQ(aux_pos[Y_AXIS])+SQ(aux_pos[Z_AXIS])-SQ(this->l2_length)-SQ(this->l3_length))/(2*this->l2_length*this->l3_length));
    ARM_Q3=atan(sqrtf(1-SQ(C3))/C3);
    S3=sinf(ARM_Q3);
    ARM_Q2=atan(aux_pos[Z_AXIS]/sqrt(SQ(aux_pos[X_AXIS])+SQ(aux_pos[Y_AXIS])))-atan((this->l3_length*sinf(ARM_Q3))/(this->l2_length+this->l3_length*cosf(ARM_Q3)));
    C2=cosf(ARM_Q2);
    S2=sinf(ARM_Q2);



    actuator_mm[ALPHA_STEPPER] = to_degrees(SCARA_theta);             // Multiply by 180/Pi  -  theta is support arm angle
    actuator_mm[BETA_STEPPER ] = to_degrees(SCARA_theta + SCARA_psi); // Morgan kinematics (dual arm)
    //actuator_mm[BETA_STEPPER ] = to_degrees(SCARA_psi);             // real scara
    actuator_mm[GAMMA_STEPPER] = cartesian_mm[Z_AXIS];                // No inverse kinematics on Z - Position to add bed offset?

}

void SixDOFArmSolution::actuator_to_cartesian(const ActuatorCoordinates &actuator_mm, float cartesian_mm[] ) const
{
    // Perform forward kinematics, and place results in cartesian_mm[]

    float y1, y2,
          actuator_rad[2];

    actuator_rad[X_AXIS] = actuator_mm[X_AXIS] / (180.0F / 3.14159265359f);
    actuator_rad[Y_AXIS] = actuator_mm[Y_AXIS] / (180.0F / 3.14159265359f);

    y1 = sinf(actuator_rad[X_AXIS]) * this->arm1_length;
    y2 = sinf(actuator_rad[Y_AXIS]) * this->arm2_length + y1;

    cartesian_mm[X_AXIS] = (((cosf(actuator_rad[X_AXIS]) * this->arm1_length) + (cosf(actuator_rad[Y_AXIS]) * this->arm2_length)) / this->morgan_scaling_x) + this->morgan_offset_x;
    cartesian_mm[Y_AXIS] = (y2 + this->morgan_offset_y) / this->morgan_scaling_y;
    cartesian_mm[Z_AXIS] = actuator_mm[Z_AXIS];

    cartesian_mm[0] = ROUND(cartesian_mm[0], 7);
    cartesian_mm[1] = ROUND(cartesian_mm[1], 7);
    cartesian_mm[2] = ROUND(cartesian_mm[2], 7);
}

bool SixDOFArmSolution::set_optional(const arm_options_t& options)
{

    arm_options_t::const_iterator i;

    i = options.find('T');         // Theta arm1 length
    if(i != options.end()) {
        arm1_length = i->second;

    }
    i = options.find('P');         // Psi arm2 length
    if(i != options.end()) {
        arm2_length = i->second;
    }
    i = options.find('X');         // Home initial position X
    if(i != options.end()) {
        morgan_offset_x = i->second;
    }
    i = options.find('Y');         // Home initial position Y
    if(i != options.end()) {
        morgan_offset_y = i->second;
    }
    i = options.find('A');         // Scaling X_AXIS
    if(i != options.end()) {
        morgan_scaling_x = i->second;
    }
    i = options.find('B');         // Scaling Y_AXIS
    if(i != options.end()) {
        morgan_scaling_y = i->second;
    }
    //i= options.find('C');          // Scaling Z_AXIS
    //if(i != options.end()) {
    //    morgan_scaling_z= i->second;
    //}
    i = options.find('D');         // Undefined min
    if(i != options.end()) {
        this->morgan_undefined_min = i->second;
    }
    i = options.find('E');         // undefined max
    if(i != options.end()) {
        this->morgan_undefined_max = i->second;
    }

    init();
    return true;
}

bool SixDOFArmSolution::get_optional(arm_options_t& options, bool force_all) const
{
    options['T'] = this->arm1_length;
    options['P'] = this->arm2_length;
    options['X'] = this->morgan_offset_x;
    options['Y'] = this->morgan_offset_y;
    options['A'] = this->morgan_scaling_x;
    options['B'] = this->morgan_scaling_y;
    // options['C']= this->morgan_scaling_z;
    options['D'] = this->morgan_undefined_min;
    options['E'] = this->morgan_undefined_max;

    return true;
};
