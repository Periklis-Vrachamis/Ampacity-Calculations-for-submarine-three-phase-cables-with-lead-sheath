-- =====================================================================
-- Submarine Three-phase Cable with Lead Sheath
-- Prepared in accordance with IEC 60287, CIGRE TB 640, and CIGRE TB 880.  
-- Verified against CIGRE TB 880, Section 6 (Case Study 2) and Section 12 (Case Study 8).
-- =====================================================================

-- =====================================================================
-- ASSUMPTIONS:
-- Conductor resistance standardized as per IEC 60228.  
-- Per CIGRE TB 880, Guidance Point 23: no lay factor is applied.  
-- T4 is valid only for cables with equal losses (see IEC 60287 §§4.2.2, 4.2.3.3.2, 4.2.3.3.3).  
-- Eddy currents are included (IEC 60287, with CIGRE TB 880 Guidance Point 6 applied).  
-- Armour losses calculated according to CIGRE TB 880, Section 2.8.  
-- =====================================================================

print("=== Submarine Three-phase Cable with Lead Sheath ===")
print("Prepared in accordance with IEC 60287, CIGRE TB 640, and CIGRE TB 880.")
print("Verified against CIGRE TB 880, Section 6 (Case Study 2) and Section 12 (Case Study 8)")

-- =====================================================================
-- SECTION 1: EXACT CABLE PARAMETERS FROM DATASHEET
-- =====================================================================

-- Cable geometry parameters (mm) 
conductor_diameter = 30                 -- Conductor diameter 
conductor_tape_diameter = 30.5          -- External diameter of conductor tape
conductor_screen_outer_diameter = 32.5  -- External diameter of sc conductor screen
insulation_outer_diameter = 48.5        -- External diameter of insulation
insulation_screen_outer_diameter = 50.3 -- External diameter of sc insulation screen
water_blocking_outer_diameter = 52.5    -- External diameter of sc water blocking tape
lead_sheath_outer_diameter = 57.1       -- External diameter of metallic (lead) sheath
pe_sheath_outer_diameter = 62.1         -- External diameter of sc PE sheath
filler_diameter = 134.9                 -- External diameter of fillers (cable after assembling)
cable_core_binder_diameter = 135.3      -- External diameter of cable core binder
armour_bedding_diameter = 137.3         -- External diameter of armour bedding
armour_outer_diameter = 149.3           -- External diameter of the armour
cable_outer_diameter = 155.3            -- External diameter of outer covering

-- Armour wires Specificatons
material_armour = "Galvanized Steel" -- Set as "Galvanized Stainless Steel" OR "Galvanized Steel"
number_armour_wires = 71 -- Number of Armour wires
diameter_armour_wires = 6  -- Diameter of Armour wires (mm)
armour_lay_length = 1785 -- Lay length of the Armour (mm)

-- System parameters 
burial_depth = 1000                    -- Burial depth to cable center (mm) 
seabed_temperature = 15.0              -- Seabed soil temperature (°C)
seabed_thermal_resistivity = 0.7       -- Soil thermal resistivity (K·m/W)
number_cables = 1                      -- Number of cables. Maximum 3 
axial_separation  = 0                  -- Axial separation between cables (mm)
max_conductor_temp = 90.0              -- Conductor operating temperature (°C)
system_voltage = 30000                -- System voltage (V)
system_frequency = 50                  -- System frequency (Hz)
R_dc_20 = 0.0283e-3                    -- Ω/m 
alpha = 3.93e-3                        -- Temperature coefficient based on IEC 60287-1-1 Table 1
skin_effect_factor = 1                 -- Based on IEC 60287-1-1 Table 2
proximity_effect_factor = 1            -- Based on IEC 60287-1-1 Table 2
epsilon = 2.5                          -- Relative permittivity of XLPE based on IEC 60287-1-1 Table 3
tan_delta = 4e-3                       -- Loss factor of XLPE (filled and unfilled) based on IEC 60287-1-1 Table 3
resistivity_sheath_20 = 21.4e-8        -- Ω·m at 20°C (electrical resistivity) based on IEC 60287-1-1 Table 1
alpha_sheath = 4.0e-3                  -- Temperature coefficient based on IEC 60287-1-1 Table 1
resistivity_armour_20 = 13.8e-8        -- Ω·m at 20°C (electrical resistivity) based on IEC 60287-1-1 Table 1
alpha_armour = 4.5e-3                  -- Temperature coefficient based on IEC 60287-1-1 Table 1
  
-- Core layup parameters - 
core_lay_length = 2152                 -- Lay length of core stranding (mm)
phase_spacing = pe_sheath_outer_diameter  -- Distance between conductor axes 

-- Calculated radius and thicknesses
r_conductor = conductor_diameter / 2
r_conductor_tape = conductor_tape_diameter / 2
r_conductor_screen = conductor_screen_outer_diameter / 2
r_insulation = insulation_outer_diameter / 2
r_insulation_screen = insulation_screen_outer_diameter / 2
r_water_blocking = water_blocking_outer_diameter / 2
r_lead_sheath = lead_sheath_outer_diameter / 2
r_pe_sheath = pe_sheath_outer_diameter / 2
r_filler = filler_diameter / 2
r_cable_core_binder = cable_core_binder_diameter / 2
r_armour_outer = armour_outer_diameter / 2
r_cable_outer = cable_outer_diameter / 2

-- =====================================================================
-- SECTION 2: MATERIAL PROPERTIES FROM IEC 60287-2-1 Table 1
-- Rationale: Thermal resistivities (rho = 1/k, where k is conductivity) 
-- =====================================================================

-- Thermal resistivities (K·m/W) 
rho_conductor = 1.0/385.0              -- Copper (assumed)
rho_conductor_tape = 6                 -- Binding tape (CIGRE TB 880 Note 2 - Table 12-11)
rho_conductor_screen = 2.5             -- SC conductor screen (IEC 60287-2-1 Table 1 and CIGRE TB 880 Table 2-3)
rho_xlpe_insulation = 3.5              -- XLPE insulation (IEC 60287-2-1 Table 1)
rho_insulation_screen = 2.5            -- SC insulation screen (IEC 60287-2-1 Table 1 and CIGRE TB 880 Table 2-3)  
rho_water_blocking = 12.0              -- Water blocking tapes (CIGRE TB 880 Table 2-3)
rho_lead_sheath = 1.0/46.7             -- Lead sheath (assumed)
rho_pe_sheath = 2.5                    -- SC PE sheath (CIGRE TB 880 Table 2-3)
rho_filler = 6.0                       -- Fillers  (CIGRE TB 880 Table 12-11)
rho_core_binder = 6.0                  -- Binding tape (CIGRE TB 880 Note 2 - Table 12-11)
rho_armour = 1.0/45                    -- Armour (assumed)
rho_outer_covering = 6.0               -- Outer covering (CIGRE TB 880 Table 2-3)
rho_seabed = seabed_thermal_resistivity

-- =====================================================================
-- SECTION 3: CORE LAYUP EFFECT - TB 880 GUIDANCE POINT 23
-- Calculate layup factor for three-core stranded arrangement
-- Rationale: Accounts for increased path length due to stranding (IEC 60287 assumes straight, TB 880 GP23 corrects).
-- =====================================================================

function calculate_core_layup_factor()
  -- Core layup factor calculation - TB 880 Table 12-2 and IEC 60287 Annex A
  local factor_Cfl = 1.29 -- Based on IEC 60287 Table A.1 for 3 number of cores
  local layup_term = pi * (factor_Cfl * pe_sheath_outer_diameter) / core_lay_length
  local f_layup = sqrt(1 + layup_term^2)
  
  print("=== CORE LAYUP EFFECT CALCULATION ===")
  print("Layup factor: " .. f_layup)
  
  return f_layup
end

layup_factor = calculate_core_layup_factor()

-- =====================================================================
-- SECTION 4: IEC 60287 THERMAL RESISTANCES 
-- Rationale: T1, T2, T3, T4 per IEC 60287-2-1 Sections 4.1 (T1 insulation), 4.2 (T2 bedding), 4.3 (T3 serving); T4 external (Section 4.2.4 for buried). 
-- =====================================================================

function calculate_Geometric_Factor_for_T2() -- CIGRE TB 880 Guidance point 45
  local x = (armour_bedding_diameter - filler_diameter) / (2 * pe_sheath_outer_diameter)
  if x <= 0 or x > 0.15 then
    print("X out of range")
  end

  if x <= 0.03 then
    -- 2π(0.00022619 + 2.11429 x − 20.4762 x^2)
    local poly = 2 * pi * (0.00022619 + 2.11429 * x - 20.4762 * x * x)
    return poly
  else
    -- 2π(0.0142108 + 1.17533 x − 4.49737 x^2 + 10.6352 x^3)
    local poly = 2 * pi * (0.0142108+ 1.17533 * x- 4.49737 * x * x + 10.6352 * x * x * x)
    return poly
  end
end

function calculate_T4() -- IEC 60287 4.2.2, 4.2.3.3.2, 4.2.3.3.3 
  if seabed_thermal_resistivity <= 0 or burial_depth <= 0 or cable_outer_diameter <= 0 then
    print("seabed_thermal_resistivity, burial_depth, and cable_outer_diameter must be positive")
  end
  local u1 = (2.0 * burial_depth) / cable_outer_diameter
  if number_cables == 1 then
    -- T4 = (1/(2π)) * ρ * ln(u1 + sqrt(u1^2 - 1))
    return (seabed_thermal_resistivity / (2.0 * pi)) * log(u1 + sqrt(u1*u1 - 1.0))

  elseif number_cables == 2 then
    -- T4 = (1/(2π)) * ρ * { ln(u1 + sqrt(u1^2 - 1)) + 0.5 * ln(1 + (2L/s1)^2) }
    return (seabed_thermal_resistivity / (2.0 * pi)) * (log(u1 + sqrt(u1*u1 - 1.0)) + 0.5 * log(1.0 + ((2.0 * burial_depth) / axial_separation)^2))

  elseif number_cables == 3 then
    -- T4 = (1/(2π)) * ρ * { ln(u1 + sqrt(u1^2 - 1)) + ln(1 + (2L/s1)^2) }
    return (seabed_thermal_resistivity / (2.0 * pi)) * (log(u1 + sqrt(u1*u1 - 1.0)) + log(1.0 + ((2.0 * burial_depth) / axial_separation)^2))
  else
    print('Error with T4')
  end
end

function calculate_thermal_resistances()
  -- T1: Conductor to metallic screen(Sheath) - IEC 60287 2-1 section 4.1.2.1 and 4.1.3
  local T1_conductor_tape = (rho_conductor_tape / (2 * pi)) * log( conductor_tape_diameter / conductor_diameter)  
  local T1_conductor_screen = (rho_conductor_screen / (2 * pi)) * log(conductor_screen_outer_diameter / conductor_tape_diameter)
  local T1_insulation = (rho_xlpe_insulation / (2 * pi)) * log(insulation_outer_diameter / conductor_screen_outer_diameter)
  local T1_insulation_screen = (rho_insulation_screen / (2 * pi)) * log(insulation_screen_outer_diameter / insulation_outer_diameter)
  local T1_water_blocking = (rho_water_blocking / (2 * pi)) * log(water_blocking_outer_diameter / insulation_screen_outer_diameter)
  
  local T1_core = T1_conductor_tape + T1_conductor_screen + T1_insulation + T1_insulation_screen + T1_water_blocking
  local T1 = T1_core / layup_factor  
  
  -- T2: Metallic screen (Sheath) to armour - IEC 60287 section 4.1.4.2
  local T2_pe_sheath = (rho_pe_sheath / (2 * pi)) * log(pe_sheath_outer_diameter / lead_sheath_outer_diameter) -- See note in the section refering to IEC 60287 2-1 section 4.1.3
  
  local G = calculate_Geometric_Factor_for_T2()
  
  local T2_fillers = (rho_filler / (6 * pi)) * G
  local T2 = (T2_pe_sheath / (3 * layup_factor)) + T2_fillers  -- T2_pe_sheath will only have the losses from one core, while the  T2_fillers should have the losses for all phases. Therefore the number of cable (3) and the effect of layup should be compenstated into the T2_pe_sheath 
  
  -- T3: Armour to cable surface - IEC 60287 section 4.1.5.1 
  local T3 = (rho_outer_covering / (2 * pi)) * log(cable_outer_diameter / armour_outer_diameter)
  
  -- T4: External thermal resistance - see above
  local T4 = calculate_T4()
  
  print("=== TB 880 THERMAL RESISTANCES ===")
  print("T1 (conductor to sheath): " .. T1 .. " K·m/W")
  print("T2 (sheath to armour): " .. T2 .. " K·m/W")
  print("T3 (armour to cable surface): " .. T3 .. " K·m/W")
  print("T4 (external): " .. T4 .. " K·m/W")
  print("Total thermal resistance: " .. T1+T2+T3+T4 .. " K·m/W")
  
  return T1, T2, T3, T4
end

-- =====================================================================
-- SECTION 5: LOSS CALCULATIONS 
-- Rationale: Losses per IEC 60287-1-1: conductor (AC resistance with skin/proximity), dielectric, sheath (approx lambda1 from TB 880 final). Adjusted for layup.
-- =====================================================================

-- Conductor AC resistance calculation 
function calculate_conductor_losses(conductor_temp)

  local factor_due_to_armour = "None"
  -- Factor related to non-magnetic and magnetic armour wires - CIGRE TB 880 Guidance Point 325 for magnetic armour wires and IEC 60287-1-1 Section 5.1.1
  if material_armour == "Galvanized Stainless Steel" then
    factor_due_to_armour = 1
  elseif material_armour == "Galvanized Steel" then
    factor_due_to_armour = 1.5
  else
    print("Armour material is wrongly set")
  end

  -- DC resistance at operating temperature - IEC 60287 1-1 section 5.1.2
  local R_dc = R_dc_20 * (1 + alpha * (conductor_temp - 20))
  
  -- Skin effect calculation - IEC 60287-1-1 Section 5.1.3
  local x_s_squared = (8 * pi * system_frequency * 1e-7 / R_dc) * skin_effect_factor
  local x_s = sqrt(x_s_squared)
  
  local y_s = "None"
  if x_s <= 2.8 then
      y_s = x_s^4 / (192 + 0.8 * x_s^4)
  elseif x_s <= 3.8 then
      y_s = -0.136 - 0.0177 * x_s + 0.0563 * x_s * x_s
  else
      y_s = 0.354 * x_s - 0.733 
  end
  
  -- Proximity effect calculation - IEC 60287-1-1 Section 5.1.5
  local x_p_squared = (8 * pi * system_frequency * 1e-7 / R_dc) * proximity_effect_factor
  local x_p = sqrt(x_p_squared)
  local dc_over_s = conductor_diameter / phase_spacing
  local y_p = (x_p^4 / (192 + 0.8 * x_p^4)) * (dc_over_s^2) * (0.312 * (dc_over_s^2) + (1.18 / ((x_p^4 / (192 + 0.8 * x_p^4)) + 0.27))) 

  -- AC resistance per meter of core - IEC 60287-1-1 Section 5.1.1
  local R_ac = R_dc * (1 + factor_due_to_armour *( y_s + y_p))
  
  -- Adjust for layup factor in case in need to be added later
  --R_ac = layup_factor * R_ac -- Attention here 
      
  return R_ac
end

-- Dielectric losses - IEC 60287 Section 5.2 and CIGRE TB 880 Guidance Point 7: Dielectric Losses
function calculate_dielectric_losses()
  local U_0 = system_voltage / sqrt(3)  -- Voltage to earth (system design)
  
  local epsilon0 = 8.854e-12 -- F/m
  -- Capacitance per meter of core -- IEC 60287 Section 5.2
  local C_core = (2.0 * pi) * epsilon0 * epsilon / log(insulation_outer_diameter / conductor_screen_outer_diameter)
  
  -- Adjust for layup factor 
  local C = layup_factor * C_core
  
  -- Angular frequency
  local omega = 2 * pi * system_frequency
  
  -- Dielectric loss per phase - IEC 60287-1-1 Section 5.2
  local W_d = omega * C * U_0^2 * tan_delta
  
  return W_d
end

-- Lead sheath losses - IEC 60287 Section 5.3
function calculate_sheath_losses(sheath_temp , R_con_ac)
  -- Lead sheath thickness calculation
  local lead_sheath_thickness = (lead_sheath_outer_diameter - water_blocking_outer_diameter) / 2

  local factor_due_to_armour = "None"
  -- Factor related to non-magnetic and magnetic armour wires - IEC 60287-1-1 Section 5.3.2 for non-magnetic armour and CIGRE TB 880 Guidance Point 38 for magnetic armour wires
  if material_armour == "Galvanized Stainless Steel" then
    factor_due_to_armour = 1
  elseif material_armour == "Galvanized Steel" then
    factor_due_to_armour = 1.5
  else
    print("Armour material is wrongly set")
  end

  -- Sheath cross-sectional area
  local A_s = pi * lead_sheath_thickness * (lead_sheath_outer_diameter - lead_sheath_thickness)  -- mm²

  -- Sheath resistance at operating temperature - IEC 60287-1-1 Section 5.3.1
  local resistivity_sheath = resistivity_sheath_20 * (1 + alpha_sheath * (sheath_temp - 20))
  local R_s = (resistivity_sheath * layup_factor) / (A_s * 1e-6)  -- Ω/m 

  -- Calculation of the reactance X per unit length of the sheath - IEC 60287 Section 5.3.2
  local mean_diam_sheath = lead_sheath_outer_diameter - lead_sheath_thickness  
  local omega = 2 * pi * system_frequency
  local X_core =  2 * omega * 1e-7 * log(2 * phase_spacing / mean_diam_sheath ) * layup_factor

  -- Calculation of lambda 1' for circulating losses of the sheath - IEC 60287-1-1 Section 5.3.2 for non-magnetic armour and CIGRE TB 880 Guidance Point 38 for magnetic armour wires
  local RS_X = R_s / X_core 
  local lambda_1_circulating = factor_due_to_armour * R_s / ( R_con_ac * (1+ RS_X * RS_X )) 

  -- Calculation of lambda 1" for eddy currents of the sheath - IEC 60287-1-1 Section 5.3.2  - This approach differs from IEC 60287-1-1 section 5.3.2. Refer to CIGRE TB 880 Guidance Point 6
  local m = omega * 1e-7 / R_s
  local lambda_0 = 3.0 * ((m*m) / (1.0 + m*m)) * ((mean_diam_sheath / (2.0*phase_spacing))^2)
  local delta_1 = (1.14 * (m^2.45) + 0.33) * ((mean_diam_sheath / (2.0 * phase_spacing)) ^ (0.92 * m + 1.66))
  local delta_2 = 0
  local b1 = sqrt( 4 * pi * omega * 1e-7 / resistivity_sheath)
  local gs = 1 + (lead_sheath_thickness/lead_sheath_outer_diameter)^1.74 * (b1 * lead_sheath_outer_diameter * 1e-3 - 1.6)
  local lambda_1_eddy = (R_s / R_con_ac) * ( gs * lambda_0 * (1.0 + delta_1 + delta_2) + ((b1 * lead_sheath_thickness)^4) / (12.0 * 1e12) )

  -- Effect of eddy currents Cf-factor - CIGRE TB 880 Guidance Point 31 - IEC 60287-1-1 Section 5.3.6
  local M = RS_X
  local N = RS_X
  local F = (4.0 * M * M * N * N + (M + N)^2) / (4.0 * (M * M + 1.0) * (N * N + 1.0))

  -- Calculation of total loss factor - IEC 60287-1-1 Section 5.3.1
  local lambda_1 = lambda_1_circulating + F * lambda_1_eddy
  
  return lambda_1, R_s, lambda_1_circulating
end

-- Calculation of Armour losses  - CIGRE TB 880 Section 
function calculate_armour_losses(armour_temp, R_con_ac, R_sheath, lambda_1_circ)
  local lambda_2 = "None"
  local mean_diam_armour = armour_outer_diameter - diameter_armour_wires
  -- Angular frequency  
  local omega = 2 * pi * system_frequency

  if material_armour == "Galvanized Stainless Steel" then
    lambda_2 = 0 -- CIGRE TB 880 Guidance Point 36 
  elseif material_armour == "Galvanized Steel" then
    -- Calculation of cross-sectional area of the armour
    local A_armour = number_armour_wires * pi * (diameter_armour_wires/2)^2 --mm²
    -- Ratio of the lenght of the wires to the length of the cable
    local layup_armour = sqrt(1+(pi * mean_diam_armour / armour_lay_length)^2)
    -- Resistance of Armour at 20oC 
    local R_AO_temp = resistivity_armour_20 / (A_armour * 1e-6)  -- Convert A_armour mm² to m²
    local k = (1.4-1.2)/(5-2) * (diameter_armour_wires - 2) + 1.2-- Guidance Point 34 and IEC 60287-1-1 Section 5.4.3.1
    local R_AO = k * R_AO_temp * layup_armour
    local R_A = R_AO * (1 + alpha_armour * (armour_temp - 20))
    -- Calculation of the value c - CIGRE TB 880 Guidance Point 33
    local C = filler_diameter / 2 - pe_sheath_outer_diameter / 2 --mm²
    -- Calculation of factor lamba_2 of the armour - CIGRE TB 880 Guidance Point 37 - IEC 60287-1-1 Section 5.4.3.5
    local factor_lamba2 = (1 - R_con_ac * lambda_1_circ / R_sheath)
    if factor_lamba2 < 0 then 
        factor_lamba2 = 0
    end
    -- Calculation of lamba_2
    lambda_2 = 1.23 * (R_A / R_con_ac) * ((2 * C * 1e-3) / (mean_diam_armour * 1e-3))^2 * (1 / (((2.77 * R_A * 1e6) / omega)^2 + 1)) * factor_lamba2 -- IEC 60287-1-1 Section 5.4.3.3.1
  else
    print("Armour material is wrongly set")
  end
  return lambda_2
end

-- Iterative calculation of current rating (ampacity)
function calculate_current_rating()
  local T1, T2, T3, T4 = calculate_thermal_resistances()
  local W_d = calculate_dielectric_losses()
  local R_ac = calculate_conductor_losses(max_conductor_temp)
  local n = 3  -- Number of phases/cores
  local delta_theta = max_conductor_temp - seabed_temperature

  -- Assumption for first iteration
  local sheath_temp = delta_theta - 5
  local armour_temp = delta_theta - 10
  local tolerance = 0.001
  local max_iter = 20

  --Set zero values
  local I = 0.0
  local iter = 0
  local lambda_1 = 0
  local lambda_2 = 0
  local R_sheath = 0
  local diff = 0

  print("=== ITERATIVE CURRENT RATING CALCULATION (AMPACITY) ===")

  repeat
    iter = iter + 1

    lambda_1, R_sheath, lambda_1_circulating = calculate_sheath_losses(sheath_temp, R_ac)
    lambda_2 = calculate_armour_losses(armour_temp, R_ac, R_sheath, lambda_1_circulating)

    local numerator = delta_theta - W_d * (0.5 * T1 + n * (T2 + T3 + T4))
    local denominator = R_ac * T1 + n * R_ac * (1 + lambda_1) * T2 + n * R_ac * (1 + lambda_1 + lambda_2) * (T3 + T4)
    local new_I = sqrt(numerator / denominator)

    local new_sheath_temp = max_conductor_temp - (new_I * new_I * R_ac + 0.5 * W_d) * T1 -- IEC 60286-1-1 Section 5.3.1
    local new_armour_temp = max_conductor_temp - ((new_I * new_I * R_ac + 0.5 * W_d) * T1 + (new_I * new_I * R_ac * (1 + lambda_1) + W_d) * n * T2) -- IEC 60286-1-1 Section 5.4.1

    -- Full per-iteration dump in compact one-line form
    print(format("Iter %d | Dphi=%.3f °C, W_d=%.6f W/m, T1=%.6f, T2=%.6f, T3=%.6f, T4=%.6f, R=%.8f Ohm/m, n=%d, lamba1=%.6f, lamba2=%.6f | I=%.3f A | sheath_used=%.3f °C", iter, delta_theta, W_d, T1, T2, T3, T4, R_ac, n, lambda_1, lambda_2, new_I, sheath_temp))

    diff = abs(new_I - I)
    I = new_I
    sheath_temp = new_sheath_temp
    armour_temp = new_armour_temp
  until diff < tolerance or iter >= max_iter
  
  print(format("Converged Current Rating (Ampacity): %.0f A", I))

  return I, R_ac, W_d,lambda_1 , lambda_2
end



-- Run the calculations
local current, R_ac, W_d, lambda_1, lambda_2 = calculate_current_rating()

-- Calculation of losses
local n = 3 -- To get the total losses in the cable they should be multiplied with n=3

-- Conductor  losses:
local conductor_losses = current * current * R_ac 

-- Dielectric losses 
local dielectric_losses = W_d 

-- Lead Sheath losses
local lead_losses = lambda_1 * current * current * R_ac 

-- Armour losses
local armour_losses = lambda_2 * current * current * R_ac * n

print("=== CALCULATION OF LOSSES PER CABLE ===") 
print(format(" Conductor Losses per core: %.2f W/m", conductor_losses))
print(format(" Dielectric Losses per core: %.2f W/m", dielectric_losses))
print(format(" Lead Sheath Losses per core: %.2f W/m", lead_losses))
print(format(" Armour Losses per cable: %.2f W/m", armour_losses))
