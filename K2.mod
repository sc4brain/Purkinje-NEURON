COMMENT

   **************************************************
   File generated by: neuroConstruct v1.6.0 
   **************************************************

   This file holds the implementation in NEURON of the Cell Mechanism:
   K2 (Type: Channel mechanism, Model: ChannelML based process)

   with parameters: 
   /channelml/@units = Physiological Units 
   /channelml/ion[1]/@name = k 
   /channelml/ion[1]/@charge = 1 
   /channelml/ion[1]/@default_erev = -85 
   /channelml/ion[1]/@role = PermeatedSubstance 
   /channelml/ion[1]/notes = The ion which actually flows through the channel 
   /channelml/ion[2]/@name = ca 
   /channelml/ion[2]/@charge = 2 
   /channelml/ion[2]/@role = ModulatingSubstance 
   /channelml/ion[2]/notes = The channel's conductance of K is dependent on this ion's concentration 
   /channelml/channel_type/@name = K2 
   /channelml/channel_type/@density = yes 
   /channelml/channel_type/status/@value = in_progress 
   /channelml/channel_type/status/comment = Getting there... 
   /channelml/channel_type/status/contributor/name = Padraig Gleeson 
   /channelml/channel_type/notes = Ca2+ dependent K current (K2). Based on Roth et al's reimplementation of original GENESIS model in NEURON 
   /channelml/channel_type/current_voltage_relation/ohmic/@ion = k 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/@default_gmax = .39 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/gate[1]/@power = 1 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/gate[1]/state/@name = m 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/gate[1]/state/@fraction = 1 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/gate[2]/@power = 2 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/gate[2]/state/@name = z 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/gate[2]/state/@fraction = 1 
   /channelml/channel_type/hh_gate[1]/@state = m 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/conc_dependence/@name = Calcium 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/conc_dependence/@ion = ca 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/conc_dependence/@variable_name = ca_conc_m 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/conc_dependence/@min_conc = 0 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/conc_dependence/@max_conc = 0.050 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/alpha/generic_equation_hh/@expr = 25 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/beta/parameterised_hh/@type = exponential 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/beta/parameterised_hh/parameter[1]/@name = A 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/beta/parameterised_hh/parameter[1]/@value = 0.075 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/beta/parameterised_hh/parameter[2]/@name = k 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/beta/parameterised_hh/parameter[2]/@value = -0.1666666667 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/beta/parameterised_hh/parameter[3]/@name = d 
   /channelml/channel_type/hh_gate[1]/transition/voltage_conc_gate/beta/parameterised_hh/parameter[3]/@value = -25 
   /channelml/channel_type/hh_gate[2]/@state = z 
   /channelml/channel_type/hh_gate[2]/transition/voltage_conc_gate/conc_dependence/@name = Calcium 
   /channelml/channel_type/hh_gate[2]/transition/voltage_conc_gate/conc_dependence/@ion = ca 
   /channelml/channel_type/hh_gate[2]/transition/voltage_conc_gate/conc_dependence/@variable_name = ca_conc_z 
   /channelml/channel_type/hh_gate[2]/transition/voltage_conc_gate/conc_dependence/@min_conc = 0 
   /channelml/channel_type/hh_gate[2]/transition/voltage_conc_gate/conc_dependence/@max_conc = 1e-8 
   /channelml/channel_type/hh_gate[2]/transition/voltage_conc_gate/alpha/generic_equation_hh/@expr = 0.2/( (ca_conc_z*1e6) *1000) 
   /channelml/channel_type/hh_gate[2]/transition/voltage_conc_gate/beta/generic_equation_hh/@expr = 1 
   /channelml/channel_type/hh_gate[2]/transition/voltage_conc_gate/tau/generic_equation_hh/@expr = 10 
   /channelml/channel_type/hh_gate[2]/transition/voltage_conc_gate/inf/generic_equation_hh/@expr = 1/(1 + alpha) 
   /channelml/channel_type/impl_prefs/comment = Using the NEURON mod file impl settings to get a closer match 
   /channelml/channel_type/impl_prefs/table_settings/@max_v = 100 
   /channelml/channel_type/impl_prefs/table_settings/@min_v = -100 
   /channelml/channel_type/impl_prefs/table_settings/@table_divisions = 200 

// File from which this was generated: /Users/hashmup/nC_projects/PurkinjeCell.ncx/cellMechanisms/K2/K2_Chan.xml

// XSL file with mapping to simulator: /Users/hashmup/nC_projects/PurkinjeCell.ncx/cellMechanisms/K2/ChannelML_v1.8.1_NEURONmod.xsl

ENDCOMMENT


?  This is a NEURON mod file generated from a ChannelML file

?  Unit system of original ChannelML file: Physiological Units

TITLE Channel: K2

COMMENT
    Ca2+ dependent K current (K2). Based on Roth et al's reimplementation of original GENESIS model in NEURON
ENDCOMMENT


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (um) = (micrometer)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (l) = (liter)
}


    
NEURON {
      

    SUFFIX K2
    USEION k READ ek WRITE ik VALENCE 1 ? reversal potential of ion is read, outgoing current is written
            
    USEION ca READ cai VALENCE 2 ? internal concentration of ion is read
            
    RANGE gmax, gion
    
    RANGE minf, mtau
    RANGE zinf, ztau
}

PARAMETER { 
      

    gmax = 0.00039 (S/cm2) ? default value, should be overwritten when conductance placed on cell
    
}



ASSIGNED {
      

    v (mV)
    
    celsius (degC)
    
    ? Reversal potential of k
    ek (mV)
    ? The outward flow of ion: k calculated by rate equations...
    ik (mA/cm2)
            
    ? The internal concentration of ion: ca is used in the rate equations...
    cai (mM)           
            
    
    gion (S/cm2)
    minf
    mtau (ms)
    zinf
    ztau (ms)
    
}

BREAKPOINT { SOLVE states METHOD derivimplicit     

    gion = gmax*((1*m)^1)*((1*z)^2)
    ik = gion*(v - ek)
                

}



INITIAL {
    ek = -85
        
    settables(v,cai)
    m = minf
        z = zinf
        
    
}
    
STATE {
    m
    z
    
}

DERIVATIVE states {
    settables(v,cai)
    m' = (minf - m)/mtau
    z' = (zinf - z)/ztau
    
}

PROCEDURE settables(v(mV), cai(mM)) {  
    
    ? Note: not all of these may be used, depending on the form of rate equations
    LOCAL  alpha, beta, tau, inf, gamma, zeta, ca_conc_m, ca_conc_z, temp_adj_m, A_alpha_m, k_alpha_m, d_alpha_m, A_beta_m, k_beta_m, d_beta_m, temp_adj_z, A_alpha_z, k_alpha_z, d_alpha_z, A_beta_z, k_beta_z, d_beta_z, A_tau_z, k_tau_z, d_tau_z, A_inf_z, k_inf_z, d_inf_z
    
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_z = 1
    
    ? Gate depends on the concentration of ca
    ca_conc_m = cai ? In NEURON, the variable for the concentration  of ca is cai
    
        
    ?      ***  Adding rate equations for gate: m  ***
             

    ? Found a generic form of the rate equation for alpha, using expression: 25
                    
    ? Equations can depend on concentration. NEURON uses 'SI Units' internally for concentration, 
    ? but ChannelML file is in Physiological Units...
    ca_conc_m = ca_conc_m / 1000000
    alpha = 25
        
    ? Resetting concentration...
    ca_conc_m = ca_conc_m * 1000000
    
    
    ? Found a parameterised form of rate equation for beta, using expression: A*exp(k*(v-d))
    A_beta_m = 0.075
    k_beta_m = -0.1666666667
    d_beta_m = -25
     
    
    beta = A_beta_m * exp((v - d_beta_m) * k_beta_m)
    
    mtau = 1/(temp_adj_m*(alpha + beta))
    minf = alpha/(alpha + beta)
          
       
    
    ?     *** Finished rate equations for gate: m ***
    
        
    ? Gate depends on the concentration of ca
    ca_conc_z = cai ? In NEURON, the variable for the concentration  of ca is cai
    
        
    ?      ***  Adding rate equations for gate: z  ***
             

    ? Found a generic form of the rate equation for alpha, using expression: 0.2/( (ca_conc_z*1e6) *1000)
                    
    ? Equations can depend on concentration. NEURON uses 'SI Units' internally for concentration, 
    ? but ChannelML file is in Physiological Units...
    ca_conc_z = ca_conc_z / 1000000
    alpha = 0.2/( (ca_conc_z*1e6) *1000)
        
    ? Resetting concentration...
    ca_conc_z = ca_conc_z * 1000000
    
         

    ? Found a generic form of the rate equation for beta, using expression: 1
                    
    ? Equations can depend on concentration. NEURON uses 'SI Units' internally for concentration, 
    ? but ChannelML file is in Physiological Units...
    ca_conc_z = ca_conc_z / 1000000
    beta = 1
        
    ? Resetting concentration...
    ca_conc_z = ca_conc_z * 1000000
    
         

    ? Found a generic form of the rate equation for tau, using expression: 10
                    
    ? Equations can depend on concentration. NEURON uses 'SI Units' internally for concentration, 
    ? but ChannelML file is in Physiological Units...
    ca_conc_z = ca_conc_z / 1000000
    tau = 10
        
    ? Resetting concentration...
    ca_conc_z = ca_conc_z * 1000000
    
    ztau = tau/temp_adj_z
         

    ? Found a generic form of the rate equation for inf, using expression: 1/(1 + alpha)
                    
    ? Equations can depend on concentration. NEURON uses 'SI Units' internally for concentration, 
    ? but ChannelML file is in Physiological Units...
    ca_conc_z = ca_conc_z / 1000000
    inf = 1/(1 + alpha)
        
    ? Resetting concentration...
    ca_conc_z = ca_conc_z * 1000000
    
    zinf = inf
          
       
    
    ?     *** Finished rate equations for gate: z ***
    
             

}


UNITSON

