#include "controller.h"
#include <iostream>
using namespace std;


controller::controller()
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        // Pos PID
        Kp_pos[i] = 0.0;
        Kd_pos[i] = 0.0;

        error_pos[i] = 0.0;
        error_old_pos[i] = error_pos[i];
        error_dot_pos[i] = 0.0;
        error_dot_old_pos[i] = error_dot_pos[i];
        PID_output_pos[i] = 0;
        PID_output_old_pos[i] = 0;

        // Vel PID
        Kp_vel[i] = 0.0;
        Kd_vel[i] = 0.0;

        error_vel[i] = 0.0;
        error_old_vel[i] = error_vel[i];
        error_dot_vel[i] = 0.0;
        error_dot_old_vel[i] = error_dot_vel[i];
        PID_output_vel[i] = 0;

        //Admittance
        deltaPos[i] = 0.0;
        deltaPos_old[i] = deltaPos[i];
        deltaPos_old2[i] = deltaPos_old[i];
        deltadot_tustin[i] = 0.0;
        deltadot_tustin_old[i] = deltadot_tustin[i];

        //RWDOB
        rhs_dob[i] = 0.0;
        rhs_dob_old[i] = rhs_dob[i];
        lhs_dob[i] = 0.0;
        lhs_dob_old[i] = lhs_dob[i];
        lhs_dob_LPF[i] = 0.0;
        lhs_dob_LPF_old[i] = lhs_dob_LPF[i];
        rhs_dob_LPF[i] = 0.0;
        rhs_dob_LPF_old[i] = rhs_dob_LPF[i];
        tauDist_hat[i] = 0.0;

        //RWFOB
        rhs_fob[i] = 0.0;
        rhs_fob_old[i] = rhs_fob[i];
        lhs_fob_LPF[i] = 0.0;
        lhs_fob_LPF_old[i] = lhs_fob_LPF[i];
        rhs_fob_LPF[i] = 0.0;
        rhs_fob_LPF_old[i] = rhs_fob_LPF[i];
        tauExt_hat[i] = 0.0;
        forceExt_hat[i] = 0.0;
        forceExt_hat_old[i] = forceExt_hat[i];
        forceExt_hat_old2[i] = forceExt_hat_old[i];

        //RWMOB
        p_hat[i] = 0.0;
        p_hat_old[i] = p_hat[i];
        p_dot_hat[i] = 0.0;
        p_dot_hat_old[i] = p_dot_hat[i];
        n_hat[i] = 0.0;
        rhs_mob[i] = 0.0;
        r[i] = 0.0;
        forceExt_MOB_hat[i] = 0.0;
        forceExt_MOB_hat_old[i] = forceExt_MOB_hat[i];
        forceExt_MOB_hat_old2[i] = forceExt_MOB_hat_old[i];
    }

};
controller::~controller(){}

void controller::pid_gain_pos(double kp, double kd, double cut_off)
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        Kp_pos[i] = kp;
        Kd_pos[i] = kd;
        cut_off_D_pos = cut_off;
    }
};

void controller::pid_gain_vel(double kp, double kd, double cut_off)
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        Kp_vel[i] = kp;
        Kd_vel[i] = kd;
        cut_off_D_vel = cut_off;
    }
};

void controller::ctrl_update()
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        //PID pos
        error_old_pos[i] = error_pos[i];
        error_dot_old_pos[i] = error_dot_pos[i];
        PID_output_old_pos[i] = PID_output_pos[i];
        //PID vel
        error_old_vel[i] = error_vel[i];
        error_dot_old_vel[i] = error_dot_vel[i];

        //admittance
        deltaPos_old2[i] = deltaPos_old[i];
        deltaPos_old[i] = deltaPos[i];
        deltadot_tustin_old[i] = deltadot_tustin[i];

        //DOB
        rhs_dob_old[i] = rhs_dob[i];
        lhs_dob_old[i] = lhs_dob[i];
        rhs_dob_LPF_old[i] = rhs_dob_LPF[i];
        lhs_dob_LPF_old[i] = lhs_dob_LPF[i];

        //FOB
        rhs_fob_old[i] = rhs_fob[i];
        rhs_fob_LPF_old[i] = rhs_fob_LPF[i];
        lhs_fob_LPF_old[i] = lhs_fob_LPF[i];
        forceExt_hat_old2[i] = forceExt_hat_old[i];
        forceExt_hat_old[i] = forceExt_hat[i];

        //MOB
        p_hat_old[i] = p_hat[i];
        p_dot_hat_old[i] = p_dot_hat[i];
        forceExt_MOB_hat_old2[i] = forceExt_MOB_hat_old[i];
        forceExt_MOB_hat_old[i] = forceExt_MOB_hat[i];
    }
};

Vector2d controller::PID_pos(StateModel_* state_model)
{
    for (int i = 0; i < NDOF_LEG; i++) // Error를 state 모델에 넣을 필요 있는지 생각해봐야함. error는 여기에 있어도 됨. //error들 update 해줘야함
    {
        error_pos[i] = state_model->posRW_ref[i] - state_model->posRW[i];
        error_old_pos[i] = state_model->posRW_ref_old[i] - state_model->posRW_old[i];

        error_dot_pos[i] = tustin_derivative(error_pos[i], error_old_pos[i], error_dot_old_pos[i], cut_off_D_pos);

        PID_output_pos[i] = Kp_pos[i] * error_pos[i] + Kd_pos[i] * error_dot_pos[i];
    }
    return PID_output_pos;
};

void controller::PID_vel(StateModel_* state_model)
{
    for (int i = 0; i < NDOF_LEG; i++)
    {
        error_vel[i] = state_model->velRW_ref[i] - state_model->velRW[i]; //error is no vector
        error_old_vel[i] = state_model->velRW_ref_old[i] - state_model->velRW_old[i];

        error_dot_vel[i] = tustin_derivative(error_vel[i], error_old_vel[i], error_dot_old_vel[i], cut_off_D_vel);

        PID_output_vel[i] = Kp_vel[i] * error_vel[i] + Kd_vel[i] * error_dot_vel[i];
    }
}; // negative velocity PID feedback

void controller::admittanceCtrl(StateModel_* state_model, double omega_n , double zeta, double k, int flag, int flag_GRF)
{
    // 현재 omega_n, zeta, k 로 tunning 하고 있는데, 변환식을 통해 아래에 적어주면 된다
    ad_M = k/(pow(omega_n,2));
    ad_B = 2*zeta*k/omega_n;
    ad_K = k;

    double c1 = 4 * ad_M + 2 * ad_B * Ts + ad_K * pow(Ts, 2);
    double c2 = -8 * ad_M + 2 * ad_K * pow(Ts, 2);
    double c3 = 4 * ad_M - 2 * ad_B * Ts + ad_K * pow(Ts, 2);

    if (flag_GRF == 0)
    {
      // using FOB
      deltaPos[0] =
          (pow(Ts, 2) * forceExt_hat[0] + 2 * pow(Ts, 2) * forceExt_hat_old[0] +
              pow(Ts, 2) * forceExt_hat_old2[0] - c2 * deltaPos_old[0] - c3 * deltaPos_old2[0]) / c1;
    }
    else if (flag_GRF == 1)
    {
      // using MOB
      deltaPos[0] =
          (pow(Ts, 2) * forceExt_MOB_hat[0] + 2 * pow(Ts, 2) * forceExt_MOB_hat_old[0] +
              pow(Ts, 2) * forceExt_MOB_hat_old2[0] - c2 * deltaPos_old[0] - c3 * deltaPos_old2[0]) / c1;
    }

    if (flag == true)
        state_model->posRW_ref[0] = state_model->posRW_ref[0] + deltaPos[0];
};

Vector2d controller::DOBRW(StateModel_* state_model, double cut_off ,int flag)
{
    cut_off_dob = cut_off;
    lhs_dob = state_model->tau_bi_old;
    rhs_dob = state_model->Lamda_nominal_DOB * state_model->qddot_bi_tustin;

    if (flag == true)
    {
        for (int i = 0; i < NDOF_LEG; i++)
        {
            lhs_dob_LPF[i] = lowpassfilter(lhs_dob[i], lhs_dob_old[i], lhs_dob_LPF_old[i], cut_off_dob);
            rhs_dob_LPF[i] = lowpassfilter(rhs_dob[i], rhs_dob_old[i], rhs_dob_LPF_old[i], cut_off_dob);

            tauDist_hat[i] = lhs_dob_LPF[i] - rhs_dob_LPF[i];
        }
    }
    else
    {
        for (int i = 0; i < NDOF_LEG; i++)
            tauDist_hat[i] = 0;
    }
    return tauDist_hat;
}; // Rotating Workspace DOB

void controller::FOBRW(StateModel_* state_model, double cut_off)
{
    cut_off_fob = cut_off;

    rhs_fob = state_model->Lamda_nominal_FOB * state_model->qddot_bi_tustin + state_model->H;

    for (int i = 0; i < NDOF_LEG; i++)
    {
        lhs_fob_LPF[i] = lowpassfilter(state_model->tau_bi_old[i], state_model->tau_bi_old2[i], lhs_fob_LPF_old[i], cut_off_fob);
        rhs_fob_LPF[i] = lowpassfilter(rhs_fob[i], rhs_fob_old[i], rhs_fob_LPF_old[i], cut_off_fob);

        tauExt_hat[i] = rhs_fob_LPF[i] - lhs_fob_LPF[i];

    }
    forceExt_hat = state_model->jacbRW_trans_inv * tauExt_hat;

    for (int i = 0; i < NDOF_LEG; i++)
    {
      state_model->GRF_FOB[i] = forceExt_hat[i];
    }
} // Rotating WorkspaceForce Observer

void controller::MOBRW(StateModel_* state_model, double ko)
{
    rhs_mob = state_model->Lamda_nominal_MOB * state_model->qdot_bi_tustin;
    n_hat = - state_model->Lamda_nominal_MOB_dot * state_model->qdot_bi_tustin + state_model->H;

    for (int i = 0; i < NDOF_LEG; i++)
    {
        p_dot_hat[i] = state_model->tau_bi_old[i] - n_hat[i] + r[i];
        p_hat[i] = tustin_integrate(p_dot_hat[i], p_dot_hat_old[i], p_hat_old[i]) + state_model->p_init[i];
        r[i] = ko * (rhs_mob[i] - p_hat[i]);
        // ko is bigger, r becomes tauExt_MOB_hat
    }
    forceExt_MOB_hat = state_model->jacbRW_trans_inv * r;

    for (int i = 0; i < NDOF_LEG; i++)
    {
      state_model->GRF_MOB[i] = forceExt_MOB_hat[i];
    }
} // Rotating WorkspaceMomentum Observer by MCL Intern Program by JihwanLee
