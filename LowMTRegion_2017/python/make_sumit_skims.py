#! /usr/bin/python
#Code to make a script to submit jobs to the localgrid
#This code doesn't submit directly jobs to localgrid, it just produces a bunch of *.sh files
#In order to submit the jobs, you'll need to type the command "source Submit_*.sh"
#Note : your code needs to be already compiler (with a .exe extension for this to work)

import os
from datasets_deepTau import * #imports dataset paths: pnfn[], myname[] and myoption[] arrays 
if __name__ == "__main__":
    location=os.getcwd();
    #name of your *compiled* code (omit the .exe extension)
    code_area = "/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion_2017/python/"
    code_name = "MuJetAnalysisBGZ.py"
    folder = "/pnfs/iihe/cms/store/user/amkalsi/Jan2019"
    cms_rel = "/user/amkalsi/CMSSW_8_0_17/src/Plots/LowMTRegion_2017/python/"
    proxy = "/user/amkalsi/x509up_u$(id -u amkalsi)"
    walltime = "15:00:00"

    for jj in range(0, len(pnfn)):    
        #Main file, which you'll use to submit the jobs
        #To submit the jobs, you'll need to type "source Submit_myname.sh" in your terminal
        submit_File = open("Submit/Submit_"+myname[jj]+".sh" , 'w')
        f=os.popen("ls -t " + pnfn[jj] + "out* | sort")   #you may want to replace "outfile" with the name of the root files found in the /pnfs folder
        dir = "dcap://maite.iihe.ac.be" +  pnfn[jj]  + "/"


        files_per_job = int(filesperjob[jj])  #Set the number of root files per job
        scr_name = ""
        ligne=0
        con_nbr = -1 
        for i in f.readlines():   #f contains the root files in the /pnfs directory
            if ligne%files_per_job==0:
                con_name = "Con"+str(con_nbr)+myname[jj]+".root"
	        con_nbr += 1   
                if (ligne != 0): 
                    outFile.write("hadd -f "+con_name+" Outout*.root\n")
#                    outFile.write("cp  "+con_name+" \t" + folder+"/Out_"+myname[jj]+"/"+con_name+"\n")
                    outFile.write("cp  "+con_name+" \t " +code_area+ "/" + "/OutputFiles/" + "\n")
                    outFile.write("rm -f *.root\n")
                    outFile.close()
                scr_name="test"+str(con_nbr)+myname[jj]
#                con_nbr += 1
                outFile = open("Jobs_to_submit/"+scr_name+".sh" , 'w')
                command1 = "source $VO_CMS_SW_DIR/cmsset_default.sh " + "\n"
                command1 = command1 + "cd "+ cms_rel + "\n" #go to the top of your CMSSW release
                command1 = command1 + "eval `scram runtime -sh` " + "\n"   #load the cms environment
                command1 = command1 + "export X509_USER_PROXY=" + proxy + "\n"  #use the proxy to write on /pnfs
                command1 = command1 + "gfal-mkdir srm://maite.iihe.ac.be:"+folder+"/Out_"+myname[jj] + "\n"  #make the ouptu sub-directory
                command1 = command1 + "export scratchdir=$TMPDIR " + "\n" #temporary directory to be used for running this job
                command1 = command1 + "cd $scratchdir" + "\n"     
	#        command1 = command1 + "mkdir Reweighting" + "\n"
        #        command1 = command1 + "cp /user/amkalsi/CMSSW_8_0_17/src/Plots/FakeRate_Corrector/Reweighting/*.root Reweighting/" + "\n" 
        #        command1 = command1 + "cp " + code_area + "/" + code_name + ".exe $scratchdir/" + "\n" #copy your code to the temporary directory
                outFile.write(command1)
                #Below: command to submit one job to the localgrid
                command3 = "qsub -q localgrid@cream02 -o "+ "out_err/"+scr_name +".stdout -e " +"out_err/"+scr_name +".stderr -l walltime="+ walltime +"    Jobs_to_submit/"+scr_name+".sh" + "\n"
                submit_File.write(command3)
            #Below, the command to execute your code, use the correct syntax for your own code, with all arguments (file in, file out, etc.) in the proper order
            #Note that i[0:-1] is just the name of the root file in /pnfs
            command2 = "dccp dcap://maite.iihe.ac.be"+   i[0:-1] + " $scratchdir/" + "\n"
            bare_fname = i[len(pnfn[jj]):-1]
#					  USE THIS ------->  MuJetAnalysisBGZ.py
            command2 = command2 + "\n" + "cmsRun  " + code_area +  code_name  + "\t" + "inputFilename="+ bare_fname +  "\t" + "outFilename=Outout" + str(ligne)+myname[jj] + ".root " + "\t" + "sample="+myoption[jj]
            command2 = command2 + " \n"+ " gfal-copy -f "+ "Outout" + str(ligne)+myname[jj] + ".root "+ " \t" + "srm://maite.iihe.ac.be:"+folder+"/Out_"+myname[jj]+"/Outout" +str(ligne)+myname[jj]+ ".root"
            command2 = command2 + "\n" + "rm -f " + bare_fname
            command2 = command2 + " \n\n\n"
            outFile.write(command2)
            ligne=ligne+1
        


        con_name = "Con"+str(con_nbr)+myname[jj]+".root"
        outFile.write("hadd -f "+con_name+" Outout*.root\n")
#        outFile.write("cp  "+con_name+" \t" + folder+"/Out_"+myname[jj]+"/"+con_name+"\n")
        outFile.write("cp  "+con_name+" \t " +code_area+ "/" + "/OutputFiles/" + "\n")
        outFile.write("rm -f *.root\n")    
        outFile.close()
        submit_File.close()
