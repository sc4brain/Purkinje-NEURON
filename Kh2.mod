COMMENT

   **************************************************
   File generated by: neuroConstruct v1.6.0 
   **************************************************

   This file holds the implementation in NEURON of the Cell Mechanism:
   Kh2 (Type: Channel mechanism, Model: Template based ChannelML file)

   with parameters: 
   /channelml/@units = Physiological Units 
   /channelml/notes = ChannelML file containing a single Channel from De Schutter and Bower 1998 
   /channelml/ion/@name = k 
   /channelml/ion/@charge = 1 
   /channelml/ion/@default_erev = -85 
   /channelml/channel_type/@name = Kh2 
   /channelml/channel_type/@density = yes 
   /channelml/channel_type/status/@value = stable 
   /channelml/channel_type/status/comment = Verified equivalence of NEURON and GENESIS mapping to orig NEURON mod impl at 0.02ms dt with current pulse 
   /channelml/channel_type/status/issue = Orig NEURON impl had to be 'fixed' to include READ ek for ion, otherwise was using the internally set value of -30! 
   /channelml/channel_type/status/contributor/name = Padraig Gleeson 
   /channelml/channel_type/notes = Anomalous rectifier current. Based on Roth et al's reimplementation of original GENESIS model in NEURON 
   /channelml/channel_type/authorList/modelAuthor[1]/name = De Schutter, E. 
   /channelml/channel_type/authorList/modelAuthor[2]/name = Bower, J. 
   /channelml/channel_type/authorList/modelTranslator[1]/name = Padraig Gleeson 
   /channelml/channel_type/authorList/modelTranslator[1]/institution = UCL 
   /channelml/channel_type/authorList/modelTranslator[1]/email = p.gleeson - at - ucl.ac.uk 
   /channelml/channel_type/authorList/modelTranslator[2]/name = Jenny Davie 
   /channelml/channel_type/authorList/modelTranslator[2]/institution = UCL 
   /channelml/channel_type/authorList/modelTranslator[2]/comment = Conversion of GENESIS model to NEURON 
   /channelml/channel_type/authorList/modelTranslator[3]/name = Arnd Roth 
   /channelml/channel_type/authorList/modelTranslator[3]/institution = UCL 
   /channelml/channel_type/authorList/modelTranslator[3]/comment = Conversion of GENESIS model to NEURON 
   /channelml/channel_type/authorList/modelTranslator[4]/name = Volker Steuber 
   /channelml/channel_type/authorList/modelTranslator[4]/institution = UCL 
   /channelml/channel_type/authorList/modelTranslator[4]/comment = Conversion of GENESIS model to NEURON 
   /channelml/channel_type/authorList/modelTranslator[5]/name = Michael Hausser 
   /channelml/channel_type/authorList/modelTranslator[5]/institution = UCL 
   /channelml/channel_type/authorList/modelTranslator[5]/comment = Conversion of GENESIS model to NEURON 
   /channelml/channel_type/publication/fullTitle = De Schutter, E., and Bower, J. M. (1994). An active membrane model of the cerebellar Purkinje cell. I. Simulation of current clamps in slice. J Neurop ... 
   /channelml/channel_type/publication/pubmedRef = http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&amp;cmd=Retrieve&amp;dopt=AbstractPlus&amp;list_uids=7512629 
   /channelml/channel_type/neuronDBref/modelName = K channels 
   /channelml/channel_type/neuronDBref/uri = http://senselab.med.yale.edu/senselab/NeuronDB/channelGene2.htm#table3 
   /channelml/channel_type/current_voltage_relation/ohmic/@ion = k 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/@default_gmax = 0.3 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/rate_adjustments/q10_settings/@q10_factor = 3 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/rate_adjustments/q10_settings/@experimental_temp = 37 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/gate/@power = 1 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/gate/state/@name = m 
   /channelml/channel_type/current_voltage_relation/ohmic/conductance/gate/state/@fraction = 1 
   /channelml/channel_type/hh_gate/@state = m 
   /channelml/channel_type/hh_gate/transition/voltage_gate/tau/generic_equation_hh/@expr = 36.8 
   /channelml/channel_type/hh_gate/transition/voltage_gate/inf/parameterised_hh/@type = sigmoid 
   /channelml/channel_type/hh_gate/transition/voltage_gate/inf/parameterised_hh/@expr = A/(1 + exp(k*(v-d))) 
   /channelml/channel_type/hh_gate/transition/voltage_gate/inf/parameterised_hh/parameter[1]/@name = A 
   /channelml/channel_type/hh_gate/transition/voltage_gate/inf/parameterised_hh/parameter[1]/@value = 0.2 
   /channelml/channel_type/hh_gate/transition/voltage_gate/inf/parameterised_hh/parameter[2]/@name = k 
   /channelml/channel_type/hh_gate/transition/voltage_gate/inf/parameterised_hh/parameter[2]/@value = 0.142857142 
   /channelml/channel_type/hh_gate/transition/voltage_gate/inf/parameterised_hh/parameter[3]/@name = d 
   /channelml/channel_type/hh_gate/transition/voltage_gate/inf/parameterised_hh/parameter[3]/@value = -82 
   /channelml/channel_type/impl_prefs/comment = Note, Using the NEURON mod file impl settings to get a closer match 
   /channelml/channel_type/impl_prefs/table_settings/@max_v = 100 
   /channelml/channel_type/impl_prefs/table_settings/@min_v = -100 
   /channelml/channel_type/impl_prefs/table_settings/@table_divisions = 200 

// File from which this was generated: /Users/hashmup/nC_projects/PurkinjeCell.ncx/cellMechanisms/Kh2/Kh2_Chan.xml

// XSL file with mapping to simulator: /Users/hashmup/nC_projects/PurkinjeCell.ncx/cellMechanisms/Kh2/ChannelML_v1.8.1_NEURONmod.xsl

ENDCOMMENT


?  This is a NEURON mod file generated from a ChannelML file

?  Unit system of original ChannelML file: Physiological Units

COMMENT
    ChannelML file containing a single Channel from De Schutter and Bower 1998
ENDCOMMENT

TITLE Channel: Kh2

COMMENT
    Anomalous rectifier current. Based on Roth et al's reimplementation of original GENESIS model in NEURON
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
      

    SUFFIX Kh2
    USEION k READ ek WRITE ik VALENCE 1 ? reversal potential of ion is read, outgoing current is written
            
    RANGE gmax, gion
    
    RANGE minf, mtau
}

PARAMETER { 
      

    gmax = 0.0003 (S/cm2) ? default value, should be overwritten when conductance placed on cell
    
}



ASSIGNED {
      

    v (mV)
    
    celsius (degC)
    
    ? Reversal potential of k
    ek (mV)
    ? The outward flow of ion: k calculated by rate equations...
    ik (mA/cm2)
            
    
    gion (S/cm2)
    minf
    mtau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
         

    gion = gmax*((1*m)^1)
    ik = gion*(v - ek)
                

}



INITIAL {
    ek = -85
        
    rates(v)
    m = minf
        
    
}
    
STATE {
    m
    
}

DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
    
}

PROCEDURE rates(v(mV)) {  
    
    ? Note: not all of these may be used, depending on the form of rate equations
    LOCAL  alpha, beta, tau, inf, gamma, zeta, temp_adj_m, A_tau_m, k_tau_m, d_tau_m, A_inf_m, k_inf_m, d_inf_m
        
    TABLE minf, mtau
 DEPEND celsius
 FROM -100 TO 100 WITH 200
    
    
    UNITSOFF
    
    ? There is a Q10 factor which will alter the tau of the gates 
                 

    temp_adj_m = 3^((celsius - 37)/10)
        
    ?      ***  Adding rate equations for gate: m  ***
             

    ? Found a generic form of the rate equation for tau, using expression: 36.8
                    tau = 36.8
        
    mtau = tau/temp_adj_m
    
    ? Found a parameterised form of rate equation for inf, using expression: A / (1 + exp(k*(v-d)))
    A_inf_m = 0.2
    k_inf_m = 0.142857142
    d_inf_m = -82
     
    
    inf = A_inf_m / (exp((v - d_inf_m) * k_inf_m) + 1)
    
    minf = inf
          
       
    
    ?     *** Finished rate equations for gate: m ***
    
             

}


UNITSON

