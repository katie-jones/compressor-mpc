/*
 * Constants and coefficients for single compressor
 */
namespace comp {
namespace coeff {
const static double J = (0.4 + 0.2070) * 0.4;
const static double tau_r = 1 / 0.5 + 1;
const static Vec<12> A((Vec<12>() << 0.000299749505193654,
                                -0.000171254191089237, 3.57321648097597e-05,
                                -9.1783572200945e-07, -0.252701086129365,
                                0.136885752773673, -0.02642368327081,
                                0.00161012740365743, 54.8046725371143,
                                -29.9550791497765, 5.27827499839098,
                                0.693826282579158).finished());

const static Vec<8> C((Vec<8>() << -0.423884232813775,
                               0.626400271518973, -0.0995040168384753,
                               0.0201535563630318, -0.490814924104294,
                               0.843580880467905, -0.423103455111209,
                               0.0386841406482887).finished());

const static Vec<8> D((Vec<8>() << -0.0083454, -0.0094965,
                               0.16826, -0.032215, -0.61199, 0.94175, -0.48522,
                               0.10369).finished());

const static double m_in_c = 0.0051;
const static Vec<2> m_rec_ss_c((Vec<2>() << 0.0047, 0.0263).finished());

const static double m_out_c = 0.017;

const static Vec<3> T_ss_c((Vec<3>() << 2.5543945754982, 47.4222669576423,
                              0.6218).finished());

const static Vec<2> SD_c((Vec<2>() << 5.55, 0.66).finished());

const static double torque_drive_c = 15000;
}
}
