
# Load Data ---------------------------------------------------------------

multiomic = list()

for(dataset in 1:6){
  ## Read Data
  switch (dataset,
          {load("Output/ihmp/proteomics_data.Rda");fname = "proteomics"},
          # {load("Output/ihmp/metabolomics_C18-neg_data.Rda");fname = "C18-neg"},
          # {load("Output/ihmp/metabolomics_C8-pos_data.Rda");fname = "C8-pos"},
          {load("Output/ihmp/metabolomics_HILIC-neg_data.Rda");fname = "HILIC-neg"},
          {load("Output/ihmp/metabolomics_HILIC-pos_data.Rda");fname = "HILIC-pos"},
          {load("Output/ihmp/mgxFunctions_data.Rda");fname = "mgxFunctions"},
          {load("Output/ihmp/metatranscriptomics_data.Rda");fname = "metatranscriptomics"},
          {load("Output/ihmp/mgxTaxonomy_data.Rda");fname = "mgxTaxonomy"}
  )
  
  ## Which samples to keep 
  load(file = "Output/ihmp/multiomic_keep.Rda")
  
  ## Get Data
  df = obj$data
  md = obj$md
  
  ##Select Key Samples
  keep  =data.frame(External.ID = keep_samples$External.ID)
  df = data.frame(External.ID = rownames(df),df)
  df = left_join(keep,df)[,-1]
  rownames(df) = keep$External.ID
  md = left_join(keep,md)
  ## View Class Labels
  table(df$Status)
  multiomic[[fname]] = df
  }


View(multiomic$proteomics)

save(multiomic,file = "Output/ihmp/ihmp_multiomic_data.Rda")
