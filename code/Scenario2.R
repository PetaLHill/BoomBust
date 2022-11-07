#simulations of boom bust in a metapopulation Scenario 2
library(slimr)
library(furrr)
library(parallel)

slim_script(
  slim_block(initialize(),
             { initializeSLiMModelType("nonWF");
               defineConstant("alpha1", 1e-9);
               defineConstant("metasz",256);
               defineConstant("subp_num",slimr_template("pNum"));
               defineConstant("r",slimr_template("rep"));
               defineConstant("mig", slimr_template("mig"));
               defineConstant("migRate",asInteger(mig/100*metasz));
               defineConstant("subp_sz", asInteger(metasz/subp_num));
               defineConstant("mev",slimr_template("mev"));
               initializeSLiMOptions(nucleotideBased=T);
               defineConstant("L",200)
               initializeAncestralNucleotides("TGCAAAATCGACAAGCCGCCCGACAGGAGTCATTCGCACCTCTGGGCCAAATTACCGACGTATGGACACAGGTTGACTTCCGGTAAGCAGGCTAGATGTTGAACCTTGGATGGCCGGAATGATCCGACAACCATGAAGACGCGTTGAGAAACGTTTTGCGTTATGTCCAGCCGATGGACCGGATCGTAAAGGGACAGTAT");
               initializeMutationTypeNuc("m1", 0.5, "f", 0.0);
               m1.convertToSubstitution = T;
               m1.mutationStackPolicy = "l";
               defineConstant("mm",matrix(c(0.0, alpha1, 0.0,0.0,alpha1,0.0,0.0,0.0,0.0,0.0,0.0,alpha1,0.0,0.0,alpha1,0.0), nrow=4, ncol=4));
               initializeGenomicElementType("g1", m1, 1.0,mm);
               initializeGenomicElement(g1, 0, L-1);
               initializeRecombinationRate(0.5);
               
             }),
  slim_block(reproduction(), {
    for (j in seqLen(subp_num+1)){
      inds = sim.subpopulations[j];
      if (inds.individualCount > 0)
      {
        parents = inds.sampleIndividuals(inds.individualCount, replace = T);
        for (i in seqLen(inds.individualCount))
        {
          mum = parents[i];
          dad = inds.sampleIndividuals(1);
          child = inds.addCrossed(mum, dad);
          child.tag = sim.tag;
          sim.tag = sim.tag + 1;
       }
      }
      else {j=j+1;
      }
      self.active = 0;}
  }),
 slim_block(1, {
   for (i in seqLen(subp_num)){
     sim.addSubpop(i,subp_sz);
     subpops = sim.subpopulations;
     f = paste("C:/INPUT_LOCATION/gl_50_" + subp_num + "_" + (i+1) + ".vcf"); 
     subpops[i]%.%genomes.readFromVCF(f, m1);
     subpops[i]%.%individuals.tag = ((subp_sz*i)+1):((i + 1) * subp_sz);
   }
   nIndividuals = sum(sim.subpopulations.individualCount);
   sim.tag = nIndividuals+1;
   sim.addSubpop(100,0);
  }),
slim_block(2,501, early(), {
    inds = sim.subpopulations.individuals;
    inds[inds.age > 0]%.%fitnessScaling = 0.0;
  }),
slim_block(2,499, first(), {
   if(sim.generation %% mev == 0){
     migrants = sample(sim.subpopulations.individuals,migRate);
     p100.takeMigrants(migrants);
     dest=p100;
    }
}),
slim_block(2,500, first(), {
    if(sim.generation %% mev == 1){
      for (id in seqLen(subp_num)){
        subpops = sim.subpopulations[sim.subpopulations.id == id];
        migrants = sample(p100.individuals, asInteger(subp_sz-subpops.individualCount));
        dest = subpops;
        dest.takeMigrants(migrants);
        }
    }
  }),
slim_block(481,500, late(), {# output from generation 498, 491, 496, 481 to generation 500 when mixing every 2,5,10 and 20 generations respectively - this example is for gene flow every 20 generations. To output a single run across last 50 generations edit this line to read slim_block(450,500, late(), {
    subpops = sim.subpopulations;
    subpops.individuals.genomes.outputVCF(paste0("C:/OUTPUT_LOCATION/file_",slimr_template("pNum"),"_mig_",slimr_template("mig"),"_rep_",slimr_template("rep"),"_mev",slimr_template("mev"),"_gen",sim.generation,".vcf"), simplifyNucleotides=T, outputNonnucleotides=F);
}),
slim_block(500,late(), {
  sim.simulationFinished();
})
)->BB

df <- expand.grid(rep=1:100,  pNum = c(2,4,8,16,32,64), mig = c(1,5,10,25,50,100), mev = c(20)) # generates a template file to run 100 replicate simulations of each combination of parameters pNum = number of subpopulations (fragmentation), mig = migration rate, mev = migration frequency. To output a single run, change rep=1:100 to rep=1
df

script_temp_df <- slim_script_render(BB, template =df)

plan(multisession(workers = 4))

results <- slim_run(script_temp_df, parallel = T, throw_error = TRUE)

