import numpy as np


def CrossSectionArea(r_tank, t_1):
    # cross sectional area is needed for allot of part this makes it so that cross section
    # can be changed everywhere if incorrect
    A = np.pi * ((r_tank ** 2) - ((r_tank - t_1) ** 2))
    return A


def MassStructure(rho, r_tank, t_1, t_2, H_tank):
    A = CrossSectionArea(r_tank, t_1)
    m_endcap = rho*(2/3)*np.pi*(r_tank**3 - (r_tank-t_2)**3)
    # mass of a single endcap simply a half of the volume of the sphere with
    # outer diameter, minus the half of the volume of the sphere with the inner diameter
    m_cylinder = rho*A*(H_tank - 2 * r_tank)  # the length of the cylinder part

    m_structure_bearing = m_endcap + m_cylinder
    m_structure = 2*m_endcap + m_cylinder
    # the bearing part is what needs to be taken into account as a compresive force
    return m_structure, m_structure_bearing


def NewDimensions(v_min, E, stress_min, pressure, p_ratio, rho, stress_uts_mat):
    L_min = 0.167  # Between these reasonable value (PLACEHOLDER)
    H_max = 4.75
    L_max = H_max
    # not really possible but ele close for which it is not possible will get picked off
    r_max = 1.293  # max radius that is possible
    f_safety = 1.25

    prop_list = dict(L_tank_list=[], r_tank_list=[], t_1_list=[], t_2_list=[], m_list=[])

    AcDif = 0.001  # accuracy of the dimension change in [m]
    for L_tank in np.arange(L_min, L_max, AcDif):

        r_over_l_fac = ((2 * stress_min) / ((np.pi ** 2) * E)) ** .5
        # written buck stress around so that you fill in stress and E and get r over l
        # if this is filled in you should get a factor that works
        r_min_column = r_over_l_fac*L_tank  # selects the minimum r possible for
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
            if 2*r_tank + L_tank < H_max:
                # else larger than maximum total height is not possible

                t_1 = (pressure*r_tank*f_safety)/stress_uts_mat

                # starting value of cylindrical region t at minimum this case due to pressure
                ShellFailure = True
                while ShellFailure:
                    stress_criticals = MaxShellBuckling(pressure, E, r_tank, t_1, L_tank,
                                                        p_ratio)
                    if stress_min < stress_criticals:
                        ShellFailure = False
                        t_1 = 1.05 * t_1  # so that there is a safety margin
                        t_2 = t_1 * 0.5
                        # (important because due to pressure there is some ratio that works
                        # for not warping the whole pressure vessel is sort of assumption)

                    else:
                        t_1 += AcDif

                v_total = (4 / 3) * np.pi * ((r_tank - t_2) ** 3) + np.pi * r_tank ** 2 * (
                        H_tank - 2 * r_tank)
                if (v_total - v_min) > 0:
                    # sometimes due to a certain t there is a decrease in volume large enough
                    # to make vtotal smaller than vmin. As thin walled approx doesnt take into
                    # account thickness and its effect on radius, it is done here


                    H_tank = L_tank + 2 * r_tank
                    m_structure, m_structure_bearing = MassStructure(rho, r_tank, t_1, t_2, H_tank)
                    prop_list["m_list"].append(m_structure)
                    prop_list["L_tank_list"].append(L_tank)
                    prop_list["r_tank_list"].append(r_tank)
                    prop_list["t_1_list"].append(t_1)
                    prop_list["t_2_list"].append(t_2)
                    # add it to dictionary or class

    min_mass_pos = prop_list["m_list"].index(min(prop_list["m_list"]))
    r_tank_new = prop_list["r_tank_list"][min_mass_pos]
    H_tank_new = prop_list["L_tank_list"][min_mass_pos] + 2 * r_tank_new
    t_1_tank_new = prop_list["t_1_list"][min_mass_pos]
    t_2_tank_new = prop_list["t_2_list"][min_mass_pos]
    # this selects for the best weighted solution
    return H_tank_new, r_tank_new, t_1_tank_new, t_2_tank_new


def StressExperienced(r_tank, t_1, t_2, rho, H_tank, acc_rocket, m_fuel, pressure, rho_fuel):

    A = CrossSectionArea(r_tank, t_1)
    # area of the cross section of the tank if only the actual material counts

    m_structure, m_structure_bearing = MassStructure(rho, r_tank, t_1, t_2, H_tank)

    v_endcap = (2 / 3) * np.pi * ((r_tank - t_2) ** 3)  # not thin walled for reducing fuel
    m_fuel_endcap = v_endcap*rho_fuel

    F_axial = (m_structure_bearing + m_fuel - m_fuel_endcap) * acc_rocket  # [N]
    F_axial_bottom = (m_structure + m_fuel) * acc_rocket
    print(m_structure, "total mass of the structure")
    print(m_structure_bearing, "total mass for stress considered")
    print(F_axial, "Axial load considered")
    print(F_axial_bottom, "Axial load total (only true for bottom)")
    print(F_axial / A, "stress compress")
    print(A, "area")

    stress_axial = (pressure * r_tank)/(2*t_1) - F_axial / A
    # f axial is removed since its in compression not tension

    print(stress_axial, "stress_axial")
    return stress_axial, F_axial


def MaxBucklingStress(r_tank, t_1, E, H_tank):
    A = CrossSectionArea(r_tank, t_1)
    I_tank = 1 / 4 * np.pi * (r_tank ** 4 - (r_tank - t_1) ** 4)
    # moment of inertia of the cross section, not thin walled assumption used

    stress_criticalb = ((np.pi ** 2) * E * I_tank) / (A * (H_tank-2*r_tank) ** 2)
    # critical column buckling stress
    print(stress_criticalb, "Critical Column buckling stress")
    return stress_criticalb


def MaxShellBuckling(pressure, E, r_tank, t_1, H_tank, p_ratio):
    Q_factor = (pressure / E) * (r_tank / t_1) ** 2

    # determine lambda value and minimum k factor
    constant_k = (12/(np.pi**4)) * (((H_tank-2*r_tank) ** 4) /
                                    ((r_tank ** 2) * (t_1 ** 2))) * (1 - p_ratio ** 2)
    # constant for k factor if k factor is equal to lambda + c/
    # the graph for this consist of 2 minimum so a plus and a minus lambda (ASSUME plus lambda)
    # ask which lambda
    lmd_half = constant_k**0.5
    # can be easily received if one calculates differentiation where a minimum is equal to 0
    k_factor = lmd_half + constant_k*(1/lmd_half)
    print(lmd_half, "lamda")
    print(k_factor, "k_factor")

    stress_criticals = (1.983 - 0.983 * np.exp(-23.14 * Q_factor)) * k_factor * \
                       ((E * np.pi ** 2) / (12 * (1 - p_ratio ** 2))) *\
                       (t_1 / (H_tank-2*r_tank)) ** 2
    # critical shell buckling (e or scientific notation)
    # oh exp can return imaginary values if not positive so be aware!!
    print(stress_criticals, "Critical sheet bucklin stress")

    return stress_criticals


def PressureStress(pressure, r_tank, t_1, t_2, stress_uts_mat):
    f_safety = 1.25
    # cylindrical part
    stress_cyl = pressure*r_tank*f_safety/t_1
    # sphere part
    stress_cap = pressure * r_tank * f_safety / (2*t_2)

    if stress_uts_mat > stress_cyl and stress_uts_mat > stress_cap:
        PressureFailure = True
    else:
        PressureFailure = False

    return PressureFailure


def PressureWeightOpt(v_min, rho, pressure, stress_uts_mat):
    r_max = 1.3  # maximum due to geometry spacecraft
    AcDif = 0.001

    H_max = 4.75
    L_max = H_max
    L_min = 0.167

    min_mass_list = dict(L_tank_list=[], r_tank_list=[], t_1_list=[], t_2_list=[], m_list=[])
    f_safety = 1.25
    for L_tank in np.arange(L_min, L_max, AcDif):
        r_forv_list = np.roots([(4 / 3) * np.pi, np.pi * L_tank, 0, -v_min])  # thin walled
        for j, i in enumerate(r_forv_list.imag):
            if i == 0 and r_forv_list.real[j] > 0:
                # this is not  nice if floating point error occurs but we
                # are assuming it doesn't occur
                # there may be a negative solution but if all values are positive one root
                # needs to be positive as well
                r_min = r_forv_list.real[j]

        for r_tank in np.arange(r_min, r_max, AcDif):
            H_tank = 2 * r_tank + L_tank
            if H_tank < H_max:
                t_1 = (pressure * r_tank * f_safety) / stress_uts_mat
                t_2 = 0.5*t_1

                v_total = (4 / 3) * np.pi * ((r_tank - t_2) ** 3) + np.pi * r_tank ** 2 * (
                        H_tank - 2 * r_tank)
                if (v_total - v_min) > 0:
                    m_structure, m_structure_bearing = MassStructure(rho, r_tank, t_1, t_2,
                                                                     H_tank)

                    min_mass_list["m_list"].append(m_structure)
                    min_mass_list["L_tank_list"].append(L_tank)
                    min_mass_list["r_tank_list"].append(r_tank)
                    min_mass_list["t_1_list"].append(t_1)
                    min_mass_list["t_2_list"].append(t_2)

    min_mass_pos = min_mass_list["m_list"].index(min(min_mass_list["m_list"]))
    r_tank_new = min_mass_list["r_tank_list"][min_mass_pos]
    H_tank_new = min_mass_list["L_tank_list"][min_mass_pos] + 2 * r_tank_new
    t_1_tank_new = min_mass_list["t_1_list"][min_mass_pos]
    t_2_tank_new = min_mass_list["t_2_list"][min_mass_pos]

    return H_tank_new, r_tank_new, t_1_tank_new, t_2_tank_new


def InputVal():
    r_tank = 0.65  # radius of the center of the fuel tank [m]
    H_tank = 0.1655 + 2 * r_tank  # length of the total tank [m] not cylindrical length
    t_1 = 0.01015  # cylindrical wall thickness [m]
    t_2 = 0.005078  # end cap thickness [m] #could calculate this with p

    # Material specific inputs
    E = 104e9  # E modulus
    p_ratio = 0.31  # poisson ratio
    rho = 4429  # density [kg/m^3]

    # other inputs
    pressure = 10000000  # [Pa]
    m_fuel = 257  # mass of the fuel in one fuel tank [kg]
    acc_rocket = 6*9.81  # highest acceleration of the rocket [m/s^2] 5.5 + 0.5
    v_min = 1.37  # minimum volume needed for fuel [#m3]
    rho_fuel = m_fuel/v_min
    stress_uts_mat = 880e6  # titanium ultimate stress

    stress_axial, F_axial = StressExperienced(r_tank, t_1, t_2, rho, H_tank, acc_rocket,
                                              m_fuel, pressure, rho_fuel)

    stress_criticalb = MaxBucklingStress(r_tank, t_1, E, H_tank)
    stress_criticals = MaxShellBuckling(pressure, E, r_tank, t_1, H_tank, p_ratio)

    PressureFailure = PressureStress(pressure, r_tank, t_1, t_2, stress_uts_mat)

    if stress_axial < 0 and \
            (PressureFailure and
             (abs(stress_axial) > stress_criticalb or abs(stress_axial) > stress_criticals)):
        H_tank_new, r_tank_new, t_1_tank_new, t_2_tank_new = \
            Iterator(PressureFailure, stress_axial, stress_criticalb, v_min, E, pressure,
                     p_ratio, rho, stress_uts_mat, acc_rocket, m_fuel, rho_fuel)

    else:
        print("")
        H_tank_new, r_tank_new, t_1_tank_new, t_2_tank_new = \
            PressureWeightOpt(v_min, rho, pressure, stress_uts_mat)

        stress_axial, F_axial = StressExperienced(r_tank_new, t_1_tank_new, t_2_tank_new, rho,
                                                  H_tank_new, acc_rocket, m_fuel, pressure,
                                                  rho_fuel)
        stress_criticalb = MaxBucklingStress(r_tank_new, t_1_tank_new, E, H_tank_new)
        stress_criticals = MaxShellBuckling(pressure, E, r_tank_new, t_1_tank_new,
                                            H_tank_new, p_ratio)
        if stress_axial < 0 and (PressureFailure and (
                abs(stress_axial) > stress_criticalb or abs(stress_axial) > stress_criticals)):
            H_tank_new, r_tank_new, t_1_tank_new, t_2_tank_new = \
                Iterator(PressureFailure, stress_axial, stress_criticalb, stress_criticals,
                         v_min, E, pressure, p_ratio, rho, stress_uts_mat, acc_rocket, m_fuel,
                         rho_fuel)

    print(H_tank_new, "total length")
    print(r_tank_new, "total radius")
    print(t_1_tank_new, "thickness cylinder")
    print(t_2_tank_new, "thickness cap")

    v_total = (4 / 3) * np.pi * ((r_tank_new - t_2) ** 3) + np.pi*r_tank_new**2 * \
              (H_tank_new - 2 * r_tank_new)

    # debugging tools
    A = CrossSectionArea(r_tank_new, t_1_tank_new)
    print(A, "crosssection")
    print(v_total, "v_total")
    print(v_min, "vmin")
    print(acc_rocket, "acceleration")


def Iterator(PressureFailure, stress_axial, stress_criticalb, stress_criticals, v_min,
             E, pressure, p_ratio, rho, stress_uts_mat, acc_rocket, m_fuel, rho_fuel):
    while PressureFailure and (
            abs(stress_axial) > stress_criticalb or abs(stress_axial) > stress_criticals):
        # check if true then thickness needs to be re assesed
        H_tank_new, r_tank_new, t_1_tank_new, t_2_tank_new = NewDimensions(v_min, E,
                                                                           abs(stress_axial),
                                                                           pressure, p_ratio,
                                                                           rho, stress_uts_mat)
        # update takes into account shell and column buckling
        # just check extra if it work or not
        stress_axial, F_axial = StressExperienced(r_tank_new, t_1_tank_new, t_2_tank_new, rho,
                                                  H_tank_new, acc_rocket, m_fuel, pressure,
                                                  rho_fuel)

        stress_criticalb = MaxBucklingStress(r_tank_new, t_1_tank_new, E, H_tank_new)
        stress_criticals = MaxShellBuckling(pressure, E, r_tank_new, t_1_tank_new, H_tank_new,
                                            p_ratio)
        PressureFailure = PressureStress(pressure, r_tank_new, t_1_tank_new, t_2_tank_new,
                                         stress_uts_mat)

    return H_tank_new, r_tank_new, t_1_tank_new, t_2_tank_new
InputVal()
