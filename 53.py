import numpy as np
# preliminary dimensions

def CrossSectionArea(r_tank, t_1):
    # cross sectional area is needed for allot of part this makes it so that cross sectionc
    # can be changed everywhere if incorrect
    A = np.pi * ((r_tank ** 2) - ((r_tank - t_1) ** 2))
    return A


def MassStructure(rho, r_tank, t_1, t_2, L_tank):
    A = CrossSectionArea(r_tank, t_1)
    m_endcap = rho*(2/3)*np.pi*(r_tank**3 - (r_tank-t_2)**3)
    # mass of a single endcap simply a half of the volume of the sphere with
    # outer diameter, minus the half of the volume of the sphere with the inner diameter
    m_cylinder = rho*A*(L_tank-2*r_tank)  # the length of the cylinder part

    m_structure_bearing = m_endcap + m_cylinder
    m_structure = 2*m_endcap + m_cylinder
    # not really mstructure but this is because this the amount of force maximum on part of the structure
    return m_structure, m_structure_bearing


def NewDimensions(v_min, E, stress_min, pressure, p_ratio, rho):
    L_min = 0.5  # Between these reasonable value (PLACEHOLDER)
    L_max = 1.2
    r_max = 1.293  # max radius that is possible
    t_1_min = 0.001  # lowest reasonable chosen value

    prop_list = dict(L_tank_list=[], r_tank_list=[], t_1_list=[], t_2_list=[], m_list=[])

    AcDif = 0.01  # accuracy of the dimension in [m]
    for L_tank in np.arange(L_min, L_max, AcDif):

        r_over_l_fac = ((2 * stress_min) / ((np.pi ** 2) * E)) ** .5
        # written buck stress around so that you fill in stress and E and get r over l
        # if this is filled in you should get a factor that works
        r_min_column = r_tank = r_over_l_fac*L_tank  # selects the minimum r possible for
        r_forv_list = np.roots([(4/3)*np.pi, np.pi*L_tank, 0, -v_min])  # thin walled
        for j, i in enumerate(r_forv_list.imag):
            if i == 0 and r_forv_list.real[j] > 0:
                # this is not really nice if floating point occurs but we dont think it does
                # there may be a negative solution but if all values are positive one root
                # needs to be positive as well
                r_min_volume = r_forv_list.real[j]

        r_min = max(r_min_column, r_min_volume) * 1.05
        # if this should fail and r_min_volume is not declared something is def not right,
        # should always have been declared. This also ensures that column buckling is done
        # the 1.05 factor is that overshoots the value so that it can slowly converge on 1.05
        for r_tank in np.arange(r_min, r_max, AcDif):
            t_1 = t_1_min
            ShellFailure = True
            while ShellFailure:
                stress_criticals = MaxShellBuckling(pressure, E, r_tank, t_1, L_tank, p_ratio)
                if stress_min < stress_criticals:
                    ShellFailure = False
                    t_1 = 1.05 * t_1  # so that it also overshoots and slowly converges to this
                    t_2 = t_1 * 0.5
                    # (important because due to pressure there is some ratio that works
                    # for not warping the whole pressure vessel is sort of assumption)

                    m_structure, m_structure_bearing = MassStructure(rho, r_tank, t_1, t_2, L_tank)
                    prop_list["m_list"].append(m_structure)
                    prop_list["L_tank_list"].append(L_tank)
                    prop_list["r_tank_list"].append(r_tank)
                    prop_list["t_1_list"].append(t_1)
                    prop_list["t_2_list"].append(t_2)
                    # add it to dictionary or class
                else:
                    t_1 += AcDif

    min_mass_pos = prop_list["m_list"].index(min(prop_list["m_list"]))
    L_tank_new = prop_list["L_tank_list"][min_mass_pos]
    r_tank_new = prop_list["r_tank_list"][min_mass_pos]
    t_1_tank_new = prop_list["t_1_list"][min_mass_pos]
    t_2_tank_new = prop_list["t_2_list"][min_mass_pos]
    # this selects for the best weighted solution
    return L_tank_new, r_tank_new, t_1_tank_new, t_2_tank_new


def StressExperienced(r_tank, t_1, t_2, rho, L_tank, acc_rocket, m_fuel, pressure):

    A = CrossSectionArea(r_tank, t_1)
    # area of the cross section of the tank if only the actual material counts

    m_structure, m_structure_bearing = MassStructure(rho, r_tank, t_1, t_2, L_tank)

    v_endcap = (2 / 3) * np.pi * ((r_tank - t_2) ** 3)  # not thin walled for reducing fuel
    v_cylinder = np.pi*r_tank**2 * (L_tank - 2*r_tank)
    m_fuel_endcap = (v_endcap/(2*v_endcap+v_cylinder))*m_fuel

    F_axial = (m_structure_bearing + m_fuel - m_fuel_endcap) * acc_rocket  # [N]
    F_axial_bottom = (m_structure + m_fuel) * acc_rocket
    print(m_structure, "total mass of the structure")
    print(m_structure_bearing, "total mass for stress considered")
    print(F_axial, "Axial load considered")
    print(F_axial_bottom, "Axial load total (only true for bottom)")
    stress_axial = (pressure * r_tank)/(2*t_1) - F_axial / A
    # f axial is removed since its in compression not tension

    print(stress_axial, "stress_axial")
    return stress_axial, F_axial


def MaxBucklingStress(r_tank, t_1, E, L_tank):
    A = CrossSectionArea(r_tank, t_1)
    I_tank = 1 / 4 * np.pi * (r_tank ** 4 - (r_tank - t_1) ** 4)
    # moment of inertia of the cross section, not thin walled assumption used

    stress_criticalb = ((np.pi ** 2) * E * I_tank) / (A * L_tank ** 2)
    # critical column buckling stress
    print(stress_criticalb, "Critical Column buckling stress")
    return stress_criticalb


def MaxShellBuckling(pressure, E, r_tank, t_1, L_tank, p_ratio):
    Q_factor = (pressure / E) * (r_tank / t_1) ** 2

    # determine lambda value and minimum k factor
    constant_k = (12/(np.pi**4)) * ((L_tank**4)/((r_tank**2) * (t_1**2))) * (1 - p_ratio ** 2)
    # constant for k factor if k factor is equal to lambda + c/
    # the graph for this consist of 2 minimum so a plus and a minus lambda (ASSUME plus lambda)
    # ask which lambda
    lmd_half = constant_k**0.5
    # can be easily received if one calculates differentiation where a minimum is equal to 0
    k_factor = lmd_half + constant_k*(1/lmd_half)
    print(lmd_half, "lamda")
    print(k_factor, "k_factor")

    stress_criticals = (1.983 - 0.983 * np.exp(-23.14 * Q_factor)) * k_factor * ((E * np.pi ** 2) / (12 * (1 - p_ratio ** 2))) * (t_1 / L_tank) ** 2
    # critical shell buckling (e or scientific notation)
    # oh exp can return imaginary values if not positive so be aware!!
    print(stress_criticals, "Critical sheet bucklin stress")

    return stress_criticals


def InputVal():
    r_tank = 0.65  # radius of the center of the fuel tank [m]
    L_tank = 0.1655 + 2 * r_tank  # length of the tank [m]
    t_1 = 0.01015  # cylindrical wall thickness [m]
    t_2 = 0.005078  # end cap thickness [m] #could calculate this with p

    # Material specific inputs
    E = 104e9  # E modulus
    p_ratio = 0.31  # poisson ratio
    rho = 4429  # density [kg/m^3]

    # other inputs
    pressure = 10000000  # [Pa]
    m_fuel = 257  # mass of the fuel in one fuel tank [kg]
    acc_rocket = 6  # highest acceleration of the rocket [m/s^2] 5.5 + 0.5
    v_min = 0.914  # minimum volume needed for fuel [#m3]

    stress_axial, F_axial = StressExperienced(r_tank, t_1, t_2, rho, L_tank, acc_rocket, m_fuel, pressure)

    stress_criticalb = MaxBucklingStress(r_tank, t_1, E, L_tank)
    stress_criticals = MaxShellBuckling(pressure, E, r_tank, t_1, L_tank, p_ratio)

    while stress_axial > stress_criticalb or stress_axial > stress_criticals:
        # check if true then thickness needs to be re assesed
        L_tank_new, r_tank_new, t_1_tank_new, t_2_tank_new = NewDimensions(v_min, E, stress_axial, pressure, p_ratio, rho)
        # update takes into account shell and column buckling
        stress_axial, F_axial = StressExperienced(r_tank_new, t_1_tank_new, t_2_tank_new, rho,
                                                  L_tank_new, acc_rocket, m_fuel, pressure)
        stress_criticalb = MaxBucklingStress(r_tank_new, t_1_new, E, L_tank_new)
        stress_criticals = MaxShellBuckling(pressure, E, r_tank_new, t_1_new, L_tank_new, p_ratio)


InputVal()
