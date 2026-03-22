import numpy as np
import config


def momentum_wheel_model(u_tau_command, ):

    # To start you get what you command
    u_tau_actual = u_tau_command

    # More advanced model

    return u_tau_actual

def translational_control(u_dv_command, ):
    u_dv_command_norm = np.linalg.norm(u_dv_command)
    # You get what you command
    if u_dv_command_norm <= config.u_cmd_max:
        u_dv_actual = u_dv_command
    else:
        u_dv_hat = u_dv_command / u_dv_command_norm
        u_dv_actual = config.u_cmd_max * u_dv_hat

    # Advanced model

    return u_dv_actual