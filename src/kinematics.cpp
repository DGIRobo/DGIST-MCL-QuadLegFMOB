#include "kinematics.h"

kinematics::kinematics(){};

kinematics::~kinematics(){};

void kinematics::state_update(StateModel_* state_model)
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        // Joint
        state_model->q_old[i] = state_model->q[i];
        state_model->q_bi_old[i] = state_model->q_bi[i]; //okay
        state_model->qdot_bi_old[i] = state_model->qdot_bi[i];
        state_model->qdot_bi_tustin_old[i] = state_model->qdot_bi_tustin[i]; //okay
        state_model->qddot_bi_tustin_old[i] = state_model->qddot_bi_tustin[i]; //okay

        // Feedback - RW Kinematics
        state_model->posRW_old[i] = state_model->posRW[i];
        state_model->velRW_old[i] = state_model->velRW[i];
        state_model->posRW_ref_old2[i] = state_model->posRW_ref_old[i];
        state_model->posRW_ref_old[i] = state_model->posRW_ref[i];
        state_model->velRW_ref_old[i] = state_model->velRW_ref[i];

        // control input
        state_model->ctrl_input_RW_old[i] = state_model->ctrl_input_RW[i];
        state_model->tau_bi_old2[i] = state_model->tau_bi_old[i];
        state_model->tau_bi_old[i] = state_model->tau_bi[i];

        // sensor noise
        state_model->Noise_old[i] = state_model->Noise[i];
    }

    state_model->H_old = state_model->H;
};


void kinematics::model_param_cal(const mjModel* m, mjData* d, StateModel_* state_model)
{
    cut_off_cal = 1/(2*pi*150);

    /* Trunk Parameters */
    m_hip = 2.5;
    m_trunk_front = 0.;
    m_trunk_rear = 0.;
    m_trunk = 4 * m_hip + m_trunk_front + m_trunk_rear;

    /* Leg Parameters */
    L = 0.25;
    d_thigh = 0.11017; // local position of CoM of thigh
    d_shank = 0.12997; // local position of CoM of shank

    m_thigh = 1.017; // mass of thigh link
    m_shank = 0.143; // mass of shank link
    m_leg = m_thigh + m_shank;
    m_total = m_trunk + 4 * m_leg;

    Izz_thigh = 0.0057;     // MoI of thigh w.r.t. CoM
    Izz_shank = 8.0318e-04; // MoI of shank w.r.t. CoM

    Jzz_thigh =
        Izz_thigh + m_thigh * pow(d_thigh, 2); // MoI of thigh w.r.t. HFE
    Jzz_shank =
        Izz_shank + m_shank * pow(d_shank, 2); // MoI of thigh w.r.t. KFE

    double M1 = Jzz_thigh + m_shank * pow(L, 2);
    double M2 = m_shank * d_shank * L * cos(state_model->q[1]);
    double M12 = Jzz_shank;

    MatInertia_bi(0,0) = M1;
    MatInertia_bi(0,1) = M12;
    MatInertia_bi(1,0) = M12;
    MatInertia_bi(1,1) = M2;

    state_model->Lamda_nominal_FOB(0,0) = M1;
    state_model->Lamda_nominal_FOB(0,1) = M12;
    state_model->Lamda_nominal_FOB(1,0) = M12;
    state_model->Lamda_nominal_FOB(1,1) = M2;

    // Added for MOB
    state_model->Lamda_nominal_MOB(0,0) = M1;
    state_model->Lamda_nominal_MOB(0,1) = M12;
    state_model->Lamda_nominal_MOB(1,0) = M12;
    state_model->Lamda_nominal_MOB(1,1) = M2;

    JzzR_thigh  = Jzz_thigh + Jzz_shank + m_shank * pow(L, 2) - 2 * m_shank * d_shank * L * cos(state_model->q[1]);
    JzzR_couple = Jzz_thigh + m_shank * pow(L, 2) - Jzz_shank;
    JzzR_shank = Jzz_thigh + Jzz_shank+ m_shank * pow(L, 2) + 2 * m_shank * d_shank * L * cos(state_model->q[1]);

    MatInertia_RW(0,0) = JzzR_thigh / (4 * pow(L, 2) * pow(sin(state_model->q[1] / 2), 2));
    MatInertia_RW(0,1) = JzzR_couple / (2 * pow(L, 2) * sin(state_model->q[1]));
    MatInertia_RW(1,0) = JzzR_couple / (2 * pow(L, 2) * sin(state_model->q[1]));
    MatInertia_RW(1,1) = JzzR_shank / (4 * pow(L, 2) * pow(cos(state_model->q[1] / 2), 2));

    Inertia_DOB(0,0) = MatInertia_RW(0,0);
    Inertia_DOB(0,1) = 0;
    Inertia_DOB(1,0) = 0;
    Inertia_DOB(1,1) = MatInertia_RW(1,1);

    double check[4] = { 0 };

    state_model->Lamda_nominal_DOB = state_model->jacbRW_trans*Inertia_DOB*state_model->jacbRW;

    //Coriolis & Gravity
    state_model->H[0] = -m_shank * d_shank * L * sin(state_model->q[1]) * pow(state_model->qdot_bi_tustin[1], 2);
         //- g * (m_thigh * d_thigh + m_shank * L) * cos(state_model->q_bi[0]);

    state_model->H[1] = m_shank * d_shank * L * sin(state_model->q[1]) * pow(state_model->qdot_bi_tustin[0], 2);
         //- g * m_shank * d_shank * cos(state_model->q_bi[1]);
}; // param_model parameter

void kinematics::sensor_measure(const mjModel* m, mjData* d, StateModel_* state_model, double cut_off)
{
    cut_off_cal = cut_off; // 150

    state_model->Noise[0] = 0.001 * gaussian_rand(0, 0.1);
    state_model->Noise[1] = 0.001 * gaussian_rand(0, 0.1);

    /*** (Serial) Joint position ***/
    state_model->q[0] = d->qpos[2] + state_model->Noise[0]; // (relative) HFE angle
    state_model->q[1] = d->qpos[3] + state_model->Noise[1]; // (relative) KFE angle

    //state_model->q[0] = d->sensordata[6]; // (relative) HFE angle
    //state_model->q[1] = d->sensordata[7]; // (relative) KFE angle

    /*** Biarticular Transformation ***/
    state_model->q_bi[0] = state_model->q[0];              // (absolute) HFE angle
    state_model->q_bi[1] = state_model->q[0] + state_model->q[1];       // (absolute) HFE angle
    //state_model->q_bi[0] = d->qpos[1];              // (absolute) HFE angle
    //state_model->q_bi[1] = d->qpos[1] + d->qpos[2]; // (absolute) KFE angle

    //state_model->q_bi[0] = d->sensordata[6];                    // (absolute) HFE angle
    //state_model->q_bi[1] = d->sensordata[6] + d->sensordata[7]; // (absolute) KFE angle

    state_model->qdot_bi[0] = d->qvel[2];
    state_model->qdot_bi[1] = d->qvel[2] + d->qvel[3];

    state_model->qddot_bi[0] = d->qacc[2];
    state_model->qddot_bi[1] = d->qacc[2] + d->qacc[3];

    for (int i = 0; i < NDOF_LEG; i++)
    {
        state_model->qdot_bi_tustin[i] =
            tustin_derivative(state_model->q_bi[i], state_model->q_bi_old[i], state_model->qdot_bi_tustin_old[i],
                cut_off_cal);
        state_model->qddot_bi_tustin[i] =
            tustin_derivative(state_model->qdot_bi_tustin[i], state_model->qdot_bi_tustin_old[i],
                state_model->qddot_bi_tustin_old[i], cut_off_cal);
    }

    // Added for MOB
    double M1_dot = 0;
    double M2_dot = - m_shank * d_shank * L * sin(state_model->q[1]) * (state_model->qdot_bi_tustin[1] - state_model->qdot_bi_tustin[0]);
    double M12_dot = 0;

    state_model->Lamda_nominal_MOB_dot(0,0) = M1_dot;
    state_model->Lamda_nominal_MOB_dot(0,1) = M12_dot;
    state_model->Lamda_nominal_MOB_dot(1,0) = M12_dot;
    state_model->Lamda_nominal_MOB_dot(1,1) = M2_dot;
};

void kinematics::jacobianRW(StateModel_* state_model)
{
    /*** Rotating Workspace ***/
    state_model->jacbRW(0,0) =  L * sin(state_model->q[1] / 2);
    state_model->jacbRW(0,1) = -L * sin(state_model->q[1] / 2);
    state_model->jacbRW(1,0) =  L * cos(state_model->q[1] / 2);
    state_model->jacbRW(1,1) =  L * cos(state_model->q[1] / 2);

    state_model->jacbRW_trans = state_model->jacbRW.transpose();

    state_model->jacbRW_trans_inv = state_model->jacbRW_trans.inverse();
};

void kinematics::fwdKinematics_cal(StateModel_* state_model)
{
    state_model->posRW[0] = 2 * L * cos((state_model->q_bi[1] - state_model->q_bi[0]) / 2); // r
    state_model->posRW[1] = (state_model->q_bi[0] + state_model->q_bi[1]) / 2;              // qr

    state_model->velRW = state_model->jacbRW*state_model->qdot_bi_tustin;

};

void kinematics::state_init(const mjModel* m, mjData* d, StateModel_* state_model)
{
    state_model->q[0] = d->qpos[2];
    state_model->q[1] = d->qpos[3];

    state_model->q_bi[0] = d->qpos[2];
    state_model->q_bi[1] = d->qpos[2] + d->qpos[3];

    // RW coordinates initialization
    state_model->r0 = 2 * L * cos((state_model->q_bi[1] - state_model->q_bi[0]) / 2);

    state_model->posRW[0] = 2 * L * cos((state_model->q_bi[1] - state_model->q_bi[0]) / 2);
    state_model->posRW[1] = (state_model->q_bi[0] + state_model->q_bi[1]) / 2;

    state_model->posRW_ref[0] = 2 * L * cos((state_model->q_bi[1] - state_model->q_bi[0]) / 2);
    state_model->posRW_ref[1] = (state_model->q_bi[0] + state_model->q_bi[1]) / 2;


    state_model->touch_sensor = 0;

    state_model->Lamda_nominal_DOB(0,0)= 0.0;
    state_model->Lamda_nominal_DOB(0,1)= 0.0;
    state_model->Lamda_nominal_DOB(1,0)= 0.0;
    state_model->Lamda_nominal_DOB(1,1)= 0.0;

    // Added for MOB
    state_model->Lamda_nominal_MOB(0,0)= 0.0;
    state_model->Lamda_nominal_MOB(0,1)= 0.0;
    state_model->Lamda_nominal_MOB(1,0)= 0.0;
    state_model->Lamda_nominal_MOB(1,1)= 0.0;

    state_model->Lamda_nominal_MOB_dot(0,0)= 0.0;
    state_model->Lamda_nominal_MOB_dot(0,1)= 0.0;
    state_model->Lamda_nominal_MOB_dot(1,0)= 0.0;
    state_model->Lamda_nominal_MOB_dot(1,1)= 0.0;

    for (int i = 0; i < NDOF_LEG; i++)
    {
        // Joint coordinates [k-1] values
        state_model->q_bi_old[i] = state_model->q_bi[i];

        state_model->qdot_bi[i] = 0.0;
        state_model->qdot_bi_tustin[i] = 0.;
        state_model->qdot_bi_tustin_old[i] = state_model->qdot_bi_tustin[i];

        state_model->qddot_bi[i] = 0.;
        state_model->qddot_bi_tustin[i] = 0.;
        state_model->qddot_bi_tustin_old[i] = state_model->qddot_bi_tustin[i];

        state_model->tau_bi[i] = 0.;
        state_model->tau_bi_old[i] = state_model->tau_bi[i];
        state_model->tau_bi_old2[i] = state_model->tau_bi_old[i];

        // RW coordinates [k-1] values
        state_model->posRW_old[i] = state_model->posRW[i];
        state_model->posRW_ref_old[i] = state_model->posRW_ref[i];
        state_model->posRW_ref_old2[i] = state_model->posRW_ref_old[i];

        state_model->velRW[i] = .0;
        state_model->velRW_old[i] = state_model->velRW[i];
        state_model->velRW_ref[i] = 0.;
        state_model->velRW_ref_old[i] = state_model->velRW_ref[i];


        state_model->ctrl_input_RW[i] = 0.;
        state_model->ctrl_input_RW_old[i] = state_model->ctrl_input_RW[i];

        // Mg Trajectory
        state_model->tau_ff[i]=0.;

        state_model->GRF_FOB[i] = 0.0;
        state_model->GRF_MOB[i] = 0.0;

        state_model->Noise[i] = 0.0;
        state_model->Noise_old[i] = state_model->Noise[i];
    }
    // Added for MOB
    state_model->p_init = state_model->Lamda_nominal_MOB * state_model->qdot_bi_tustin;
};

double kinematics::gaussian_rand(double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}
