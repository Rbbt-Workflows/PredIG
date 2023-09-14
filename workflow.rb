require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/util/R'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/PredIG'

Workflow.require_workflow "Immunomics"
Workflow.require_workflow "NOAH"
Workflow.require_workflow "NetCleave"

module PredIG
  extend Workflow



  input :peptides, :array, "Array of peptides", nil, :required => true
  task :tapmat_pred_fsa => :tsv do |peptides|

    peptides = Open.read(peptides).split("\n") if String ===  peptides

    translations = {}
    fasta = ""
    peptides.each do |peptide|
      fasta += <<-EOF
>#{peptide}
#{peptide}
      EOF
    end

    sizes = peptides.collect{|p| p.length }.uniq.sort

    Open.write file('fasta'), fasta
    
    scores = TSV.setup({}, "Peptide~TAP#:type=:single#:cast=:to_f")
    sizes.each do |size|
      scores.merge!(Immunomics.tapmat_pred_fsa(file('fasta'), size))
    end

    scores.subset(peptides)
  end

  input :peptides, :array, "List of peptides, one per line, no header", nil, :required => true
  task :physico_chemical => :tsv do |peptides|

    peptides = Open.read(peptides).split("\n") if String ===  peptides
    
    R.run <<-EOF, :log => true
rbbt.require("Peptides")
rbbt.require("dplyr")
rbbt.require("xgboost")
rbbt.require("argparser")
rbbt.require("stringr")
rbbt.require("seqinr")
    EOF

    str = 'peptide' << "\n"
    str += peptides * "\n"

    input = file('input.csv')

    Open.write(input, str + "\n")
    CMD.cmd_log("Rscript #{Rbbt.share.PredIG.R["predig_pch_calc.R"].find} --input '#{input}'")

    TSV.open(input.sub('.csv','_pch.csv'), :sep => ',', :header_hash => '', :type => :list)
  end

  # ToDo: Accept the peptide and allele columns on other order
  input :input, :text, "Two columns with peptide and allele"
  task :mhcFlurry => :tsv do |input|
    str = 'allele,peptide' << "\n"

    input = StringIO.new input unless Misc.is_filename?(input)
    TSV.traverse input, :into => str, :type => :array do |line|
      next if line.downcase.include?("peptide")
      peptide, allele = line.split(/,|\t/)
      [allele, peptide]  * "," + "\n"
    end

    input = file('input.csv')
    output = file('output.csv')

    Open.write(input, str)
    CMD.cmd_log('mhcflurry-predict'," #{ input } --out #{output}" )

    mhcflurry_scores = TSV.setup({}, :key_field => "ID", :fields => ["mhcflurry_affinity", "mhcflurry_affinity_percentile","mhcflurry_processing_score", "mhcflurry_presentation_score"], :type => :list)
    TSV.traverse output, :type => :array do |line|
      next if line =~ /peptide/
      allele, peptide, pred, percent, processing, presentation = line.split(",")
      mhcflurry_scores[[peptide, allele] * "_"] = [pred, percent, processing, presentation]
    end

    mhcflurry_scores
  end

  input :input, :text, "Three columns with Peptide, UniProt_ID, and Allele"
  task :prepare_input_files => :array do |input|
    peptides = []
    peptide_alleles = []
    peptide_uniprot = []
    TSV.traverse StringIO.new(input), :type => :array do |line|
      next if line.downcase.include? 'peptide'
      peptide, allele, uniprot = line.split(/,|\t/)
      peptides << peptide
      peptide_alleles << [peptide, allele] * ","
      peptide_uniprot << [peptide, uniprot] * ","
    end

    Open.write(file('peptides'), peptides.uniq * "\n")
    Open.write(file('peptide_alleles'), "peptide,allele\n" + peptide_alleles.uniq * "\n")
    Open.write(file('peptide_uniprot'), "peptide,uniprot\n" + peptide_uniprot.uniq * "\n")
    file.glob("*")
  end

  dep :prepare_input_files
  dep :tapmat_pred_fsa, :peptides => :placeholder do |jobname,options,dependencies|
    input_file_job = dependencies.flatten.first
    {:inputs => options.merge(:peptides => input_file_job.file('peptides'))}
  end
  dep :physico_chemical, :peptides => :placeholder do |jobname,options,dependencies|
    input_file_job = dependencies.flatten.first
    {:inputs => options.merge(:peptides => input_file_job.file('peptides'))}
  end
  dep :mhcFlurry, :input => :placeholder do |jobname,options,dependencies|
    input_file_job = dependencies.flatten.first
    {:inputs => options.merge(:input => input_file_job.file('peptide_alleles'))}
  end
  dep NOAH, :predict, :input => :placeholder do |jobname,options,dependencies|
    input_file_job = dependencies.flatten.first
    {:inputs => options.merge(:input => input_file_job.file('peptide_alleles'))}
  end
  dep NetCleave, :predict_uniprot, :uniprot_csv => :placeholder do |jobname,options,dependencies|
    input_file_job = dependencies.flatten.first
    {:inputs => options.merge(:uniprot_csv => input_file_job.file('peptide_uniprot'))}
  end
  task :predig => :tsv do
    tap, pch, mhcFlurry, noah, net_cleave = dependencies[1..-1].collect{|dep| dep.load }

    noah.add_field "Peptide_Allele" do |k|
      k.split("_").reverse * "_"
    end

    noah = noah.reorder "Peptide_Allele", noah.fields - ["Peptide_Allele"]
    input = noah.attach tap
    pch.key_field = "Peptide"
    input = input.attach pch
    net_cleave = net_cleave.reorder "epitope", "netcleave"
    net_cleave.key_field = "Peptide"
    input = input.attach net_cleave
    mhcFlurry.key_field = "Peptide_Allele"
    input = input.attach mhcFlurry

    input = input.slice input.fields - %w(Allele Peptide tcr_contact)
    input.fields = input.fields.collect{|f| f == "NOAH_score" ? "NOAH" : f }
    input

    input_csv = file('input.csv')

    peptide_alleles = []
    str = input.fields * "," + "\n"
    TSV.traverse input do |k,values|
      peptide_alleles << k
      values = values.collect{|v| v.to_s}
      str += values*"," + "\n"
    end

    Open.write(input_csv, str)

    CMD.cmd_log("Rscript '#{Rbbt.share.PredIG.R["predig_spwindep_calc.R"].find}' --input '#{input_csv}' --model '#{Rbbt.share.PredIG.model["predig_imbal_grid_spw_indep.model"].find}'")

    scores = []
    Open.read(input_csv.sub('.csv', '_predig.csv')).split("\n").each do |line|
      next if line.downcase.include? "noah"
      predig_score = line.split(",").last
      scores << predig_score
    end

    score_tsv = TSV.setup(Misc.zip2hash(peptide_alleles, scores), :key_field => "Peptide_Allele", :fields => ["PredIG"], :type => :single)

    output = input.attach score_tsv

    output.add_field "Peptide" do |k|
      k.split("_").first
    end

    output.add_field "Allele" do |k|
      k.split("_").last
    end

    output.reorder :key, (["Peptide", "Allele", "PredIG"] + output.fields).uniq
  end



end

#require 'PredIG/tasks/basic.rb'

#require 'rbbt/knowledge_base/PredIG'
#require 'rbbt/entity/PredIG'

