include("Include.jl")

# set some constants -
AV_NUMBER=6.02e23



function load_parameter_file(path_to_parameter_file::String)::Dict{String,Any}
    return TOML.parsefile(path_to_parameter_file);
end

function load_experimental_data_file(path_to_data_file::String, parameter_dictionary::Dict{String,Any})::Array{Float64,2}
    
    # raw -
    raw_data = readdlm(path_to_data_file,',');
    (NR,NC) = size(raw_data)
    converted_data_array = zeros(NR,NC)

    # need to convert to nmol/gDW -
    water_fraction = parameter_dictionary["fraction_of_water_ecoli_cell"]["value"]      # units: dimensionless
    CW = parameter_dictionary["mass_of_single_ecoli_cell"]["value"]                     # units: g/cell
    gDW = (1-water_fraction)*CW                                                         # units: gDW/cell

    # scale factor -
    SF = (1/AV_NUMBER)*(1/gDW)*1e9

    # convert -
    for row_index = 1:NR
        converted_data_array[row_index,1] = raw_data[row_index,1]       # 1st col: I concentration -
        converted_data_array[row_index,2] = SF*raw_data[row_index,2]    # 2nd col:  mRNA
        converted_data_array[row_index,3] = SF*raw_data[row_index,2]    # 3rd col:  mRNA low
        converted_data_array[row_index,4] = SF*raw_data[row_index,2]    # 4th col:  mRNA high
    end
    
    # return -
    return converted_data_array 
end

function simulate_lacZ_expression(GM,W1,W2,Kd,n)

    # initialize -
    lacZ_sim_array = zeros(1000,3)

    # setup inducer range -
    I_range = exp10.(range(-5,0,length=1000))

    # compute -
    for (index,i_val) in enumerate(I_range)
        
        # grab inducer level -
        lacZ_sim_array[index,1] = i_val
        
        # compute u -
        f_val = (i_val^n)/(Kd^n+i_val^2)
        u_value = (W1+W2*f_val)/(1+W1+W2*f_val)

        # compute mRNA 
        lacZ_sim_array[index,2] = GM*u_value

        # grab u -
        lacZ_sim_array[index,3] = u_value
    end

    # return -
    return lacZ_sim_array
end

function compute_gain_function_model(parameter_dictionary::Dict{String,Any})

    # the gain is the transcriptional kinetic limit divided by the sink terms (dilution + degradation)
    
    # general stuff -
    volume_of_ecoli_cell = parameter_dictionary["volume_of_ecoli_cell"]["value"]        # units: mum^3
    volume_of_ecoli_cell_L = volume_of_ecoli_cell*((1/1e6)^3)*((100)^3)*(1/1000)        # units: L
    water_fraction = parameter_dictionary["fraction_of_water_ecoli_cell"]["value"]      # units: dimensionless
    CW = parameter_dictionary["mass_of_single_ecoli_cell"]["value"]                     # units: g/cell
    gDW = (1-water_fraction)*CW                                                         # units: gDW/cell
    frac_active_RNAP = parameter_dictionary["fraction_RNAP_active_ecoli_cell"]["value"] # units: frac active

    # ---- denominator - the degradation + dilution terms -------------------------------------------- #
    cell_double_time = parameter_dictionary["ecoli_doubling_time"]["value"]             # units: min
    mRNA_half_life = parameter_dictionary["mRNA_half_life"]["value"]                    # units: min
    specific_growth_rate = (log(2)/cell_double_time)*60.0                               # units: hr^-1
    mRNA_degradation_rate_constant = -(log(0.5)/(mRNA_half_life))*60.0                    # units: hr^-1
    denominator = (specific_growth_rate + mRNA_degradation_rate_constant)               # units: hr^-1
    # ------------------------------------------------------------------------------------------------ #

    # ---- numerator - the transcriptional limit ----------------------------------------------------- #
    # Kx -
    m = parameter_dictionary["slope_initiation_exepriment"]["value"]                    # units: muM/s
    yi = parameter_dictionary["intercept_initiation_experiment"]["value"]               # units: sec
    kI = 1/(yi)                                                                         # units: sec^-1
    KX = m*kI*(volume_of_ecoli_cell_L/gDW)*(1e9/1e6)                                    # units: nmol/gDW       
    
    # RANP concentration -
    RNAP_copy_number = parameter_dictionary["average_RNAP_copy_number_ecoli_cell"]["value"]    # units: #/cell
    RNAP_concentration = (RNAP_copy_number)*(1/AV_NUMBER)*(1/gDW)*1e9                          # units: nmol/gDW
    RNAP_active = (frac_active_RNAP)*RNAP_concentration

    # lacZ gene concentration -
    lacZ_copy_number = parameter_dictionary["lacZ_copy_number"]["value"]                # units: #/cell
    G = (lacZ_copy_number)*(1/AV_NUMBER)*(1/gDW)*1e9                                    # units: nmol/gDW
    
    # time constant -
    RNAP_elongation_rate_ecoli_cell = parameter_dictionary["RNAP_elongation_rate_ecoli_cell"]["value"]  # units: nt/s
    L = parameter_dictionary["gene_charateristic_length"]["value"]                      # units: nt
    lacZ_L = parameter_dictionary["length_lacZ_gene"]["value"]                          # units: nt
    avg_eln_constant = (RNAP_elongation_rate_ecoli_cell)*(1/L)                          # units: s^-1
    lacZ_eln_constant = avg_eln_constant*(L/lacZ_L)                                     # units: s^-1
    tau = (lacZ_eln_constant)/kI                                                        # units: dimensionless
    numerator = lacZ_eln_constant*(RNAP_active)*(G/(KX*tau+G*(1+tau)))*(3600)           # units: nmol/gDW-hr
    # ------------------------------------------------------------------------------------------------ #

    # return -
    return ((numerator/denominator), KX, tau)
end

# setup -
path_to_data_file="./data/Data.dat"
path_to_parameter_file="./Parameters.toml"

# load the parameters -
pd = load_parameter_file(path_to_parameter_file)

# load and convert the data -
data_array = load_experimental_data_file(path_to_data_file, pd);

# compute the gain model -
(GM, KX, tau) = compute_gain_function_model(pd)

# estimate W1 -
alpha = (data_array[1,2]/GM)    # m(I=0)/K = W1/(1+W1)
W1 = (alpha)/(1-alpha)

# estimate W2 -
# saturates at I = 0.216 => assume fI ~ 1 => we can solve for W2 directly u = 0.99
u_sat = 0.99
W2 = (1/(1-u_sat))*(u_sat+(u_sat - 1)*W1)

# how do we get KD and n?
# we can estimate these by least-squares (correct way), or just fiddle w/them based on estimates -
Kd = 9e-2   # units: mM
n = 1.85    # units: dimensionless
sim_array_lacZ = simulate_lacZ_expression(GM,W1,W2,Kd,n)

# make mRNA plot -
figure()
semilogx(sim_array_lacZ[:,1],sim_array_lacZ[:,2],"b",lw=2)
semilogx(data_array[:,1], data_array[:,2], "o",mec="b", mfc="white")
xlabel("Extracellular IPTG [mM]", fontsize=16)
ylabel("lacZ mRNA concentration [nmol/gDW]", fontsize=14)

# set the face color -
fig = gcf()
ax = gca()
ax.set_facecolor("#F5FCFF")
fig.set_facecolor("#ffffff")

# Make the u-plot -
figure()
semilogx(sim_array_lacZ[:,1],sim_array_lacZ[:,3],"r",lw=2)
xlabel("Extracellular IPTG [mM]", fontsize=16)
ylabel(L"Promoter model $\bar{u}$ (AU)", fontsize=14)

# write -
# savefig("./lacZ-mRNA-model-v-data.pdf", facecolor=fig.get_facecolor(), edgecolor="none")
