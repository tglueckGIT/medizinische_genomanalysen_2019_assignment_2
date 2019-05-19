#! /usr/bin/env python3

import vcf

__author__ = 'Glueck Tobias'


class Assignment2:
    
    def __init__(self):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)
        

    def get_average_quality_of_file(self):
        with open("chr22.vcf") as f:
            c = 0
            sum = 0
            for line in f:
                c += 1
                if len(line.split("\t"))>4 and line.split("\t")[5].isdigit() == True:
                    sum += int(line.split("\t")[5])
            return(sum/c)
        
        
    def get_total_number_of_variants_of_file(self):
        with open("chr22.vcf") as f:
            c = 0
            for line in f:
                if line.split("\t")[0] == "chr22":
                    c += 1
            return(c)
    
    
    def get_variant_caller_of_vcf(self):
        with open("chr22.vcf") as f:
            for line in f:
                if line.split("\t")[0] == "chr22":
                    start = line.find("callsetnames=") + 13
                    end = line.find("datasetsmissingcall=") - 1
                    return(line[start:end])

        
        
        
    def get_human_reference_version(self):
        with open("chr22.vcf") as f:
            for line in f:
                if line.split("\t")[0] == "chr22":
                    start = line.find("difficultregion=") + 16
                    if start == 15:
                        continue
                    end = line.find(";" , start) - 1
                    return(line[start:end])
        
        
    def get_number_of_indels(self):
        indels = 0

        for i in vcf.Reader(open("chr22.vcf", "r")):
            if i.is_indel:
                indels += 1
        return indels
        

    def get_number_of_snvs(self):
        snvs = 0

        for i in vcf.Reader(open("chr22.vcf", "r")):
            if i.is_snp:
                snvs += 1
        return snvs
        
    def get_number_of_heterozygous_variants(self):
        hetvar = 0

        for i in vcf.Reader(open("chr22.vcf", "r")):
                hetvar += i.num_het
        return hetvar
        
    
    def merge_chrs_into_one_vcf(self):
        f1 = vcf.Reader(open("chr21.vcf"), "r")
        f2 = vcf.Reader(open("chr22.vcf"), "r")
        merge = vcf.Writer(open("chr21_chr22_merged.vcf", "w"), f1)

        for file in [f1, f2]:
            for line in file:
                merge.write_record(line)

        with open("chr21_chr22_merged.vcf") as f:
            c = 0
            for line in f:
                if line.split("\t")[0] == "chr22":
                    c += 1
            return(c)
    
    def print_summary(self):
        print("Average quality: " + str(self.get_average_quality_of_file()))
        print("Total number of variants: " + str(self.get_total_number_of_variants_of_file()))
        print("Variant caller: " + str(self.get_variant_caller_of_vcf()))
        print("Human reference version: " + (str(self.get_human_reference_version())))
        print("Number of indels: " + str(self.get_number_of_indels()))
        print("Number of snvs: " + str(self.get_number_of_snvs()))
        print("Number of heterozygous variants: " + str(self.get_number_of_heterozygous_variants()))
        print("Total number of variants in merge file: " + str(self.merge_chrs_into_one_vcf()))
        
    
    
def main():
    print("Assignment 2")
    assignment2 = Assignment2()
    assignment2.print_summary()
    print("Done with assignment 2")
        
        
if __name__ == '__main__':
    main()
   
    



