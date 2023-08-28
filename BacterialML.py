import pandas as pd 
import keras 
import keras_dna
import tensorflow as tf
from keras_dna import Generator
from keras_dna import ModelWrapper, Generator
from tensorflow.keras.models import Sequential
class BacterialML:
    """summary_line
    a keras dna implementation of the machine learning class
    for the sequences from the bacterial genomes from the newly 
    sequenced genomes.  
    Keyword arguments:
    argument -- description
    gff_: annotation_gff of the bacterial genomes
    batch_size : batch_size
    gff_out_ : path for the gff converted files.
    genes_ : genes names on which model to train
    Return: return_description
    """
    def __init__(self, read_gff, write_gff_out, batch_size): 
        self.gff = read_gff
        self.gff_out = write_gff_out
        self.batch_size = batch_size
        with open(self.gff, "r") as read_file:
            with open(self.gff_out,"w") as write_file:
                for line in read_file.readlines():
                    if line.startswith("#"):
                        continue
                    write_file.write(line)
                self.gff.close()
                self.gff_out.close()
    def readModelSelection(self,genes):
        self.gff_model_selection = pd.read_csv("self.gff_out", sep = "\t")
        self.genes = genes
        self.genes_read = []
        with open(self.genes, "r") as gene_read:
            for line in gene_read.readlines():
                self.genes_read.append(line.strip())
                gene_read.close()
    
    def trainingModel(self, fasta, size):
        self.fasta = fasta
        self.size = size
        self.size_estimate = list(filter(None,[x.strip() for x in open(self.fasta).readlines()])) 
        self.size_estimates = []
        fasta_read = {}
        for i in self.size_estimate:   
            if i.startswith(">"):
                fasta_read_path = i.strip()
            if i not in fasta_read:
                fasta_read[i] = ""
                continue
            fasta_read[fasta_read_path] += i.strip()
        fasta_ids = list(fasta_read.keys())
        fasta_sequences =  list(fasta_read.values())
        fasta_dataframe = pd.DataFrame([(i,j)for i,j in zip(fasta_ids, fasta_sequences)]). \
                                          rename(columns = {0: "ids", 1: "sequence"})
                                          
        fasta_dataframe["length"] = fasta_dataframe["sequences"].apply(lambda n: len(n))
        self.size_estimates = fasta_dataframe["length"].to_list()                               
        bacterial_model_generator = Generator(batch_size = int(self.batch_size),
                                           fasta_file = "self.fasta",
                                           annotation_files = ["self.gff_out"],
                                           annotation_list = [self.genes_read[i] \
                                                              for i in range(len(self.genes_read))])
        bacterial_model = Sequential()
        bacterial_model.compile(loss = "mse", optimizer = "adam")
        bacterial_wrap = ModelWrapper(model = bacterial_model, \
                                       generator_train = bacterial_model_generator)
        bacterial_model.wrap.train(epochs = 10)
        bacterial_model.evaluate(incl_chromosome = ["chr"])
        bacterial_model.predict(incl_chromosome = ["chr"], chrom_size = "self.size")
        # if you want to take all the chromosome sizes then 
        bacterial_model.predict(incl_chromosome = ["chr"], chrom_size = "self.size.estimate")
        bacterial_model.wrap.save(path = "self.__path", save_model = True)
