def Condor(runscript):
    dir_base = os.environ['CMSSW_BASE']+'/src/2DAlphabet/'

    commands = []

    # tar_cmd = "tar czvf tarball.tgz --directory="+dir_base+" "
    # for t in tarFilesList: tar_cmd+=t+' '
    # commands.append(tar_cmd)

    # Make JDL from template
    timestr = time.strftime("%Y%m%d-%H%M%S")
    out_jdl = 'temp_'+timestr+'_jdl'
    commands.append("sed 's$TEMPSCRIPT$"+runscript+"$g' "+dir_base+"condor/templates/jdl_template_noargs > "+out_jdl)
    commands.append('mkdir notneeded')
    commands.append("condor_submit "+out_jdl)
    commands.append("mv "+out_jdl+" "+dir_base+"condor/jdls/")

    # commands.append("condor_q lcorcodi")

    for s in commands:
        print s
        subprocess.call([s],shell=True)

def StatsForCondor(run_name,toyDict,tarFilesList,commands,files_to_grab=[]):
    if len(files_to_grab) == 0:
        files_to_grab = ['higgsCombine'+run_name+'*.*.mH120.*.root','fitDiagnostics'+run_name+'*.root']
    ntoys = toyDict['ntoys']
    toys_per_job = toyDict['toys_per_job']
    jobs = toyDict['toyjobs']
    seed = toyDict['seed']

    executeCmd('mkdir condor_'+run_name)
    executeCmd('rm condor_'+run_name+'/*')

    with cd('condor_'+run_name):
        # Ship along a post-processing script to grab the output from condor
        tar_cmd = "tar czvf tarball.tgz --directory="+os.environ['CMSSW_BASE']+'/src/2DAlphabet/'+" "
        for t in tarFilesList: tar_cmd+=t+' '
        executeCmd(tar_cmd)

        out_tar_cmd = 'tar -czvf '+run_name+'.tgz '
        for f in files_to_grab: out_tar_cmd+=f+' '

        if jobs > 1:
            for i in range(jobs):
                executeCmd('cp '+os.environ['CMSSW_BASE']+'/src/2DAlphabet/condor/templates/run_combine.sh run_combine_%s.sh'%(i))
                this_run_file = open('run_combine_%s.sh'%(i),'a+')

                this_seed = random.randint(100000,999999)
                this_run_name = run_name+'_%s'%(this_seed)

                these_commands = []
                for cmd in commands:
                    this_command = cmd.replace('-t '+str(ntoys),'-t '+str(toys_per_job))
                    this_command = this_command.replace('-n '+run_name,'-n '+this_run_name)
                    this_command = this_command.replace('-s '+str(seed),'-s '+str(this_seed))
                    # this_command = this_command.replace('-d '+workspace_name,'-d '+projDir+'/'+workspace_name)
                    this_command = this_command.replace('.'+str(seed)+'.','.'+str(this_seed)+'.')
                    this_command = this_command.replace('higgsCombine'+run_name,'higgsCombine'+this_run_name)
                    these_commands.append(this_command)

                this_run_file.write('\n')
                for c in these_commands:
                    this_run_file.write(c+'\n')
                this_run_file.write(out_tar_cmd.replace(run_name+'.tgz',run_name+'_'+str(this_seed)+'.tgz')+'\n')
                this_run_file.write('cp '+this_run_name+'.tgz $CMSSW_BASE/../')
                this_run_file.close()
                Condor('run_combine_%s.sh'%(i))
                
        else:
            executeCmd('cp '+os.environ['CMSSW_BASE']+'/src/2DAlphabet/condor/templates/run_combine.sh run_combine.sh')
            this_run_file = open('run_combine.sh','a+')
            this_run_file.write('\n')
            for c in commands:
                this_run_file.write(c+'\n')
            this_run_file.write(out_tar_cmd+'\n')
            this_run_file.write('cp '+run_name+'.tgz $CMSSW_BASE/../')
            this_run_file.close()
            Condor('run_combine.sh')

def MakeShellFinisher(identifier,files_to_grab):
    tar_cmd = 'tar -czvf '+identifier+'.tgz '
    for f in files_to_grab: tar_cmd+=f+' '
    shell_finisher = open('shell_finisher.sh','w')
    shell_finisher.write('#!/bin/sh\n')
    shell_finisher.write('echo "'+tar_cmd+'" \n')
    shell_finisher.write(tar_cmd+'\n')
    shell_finisher.write('cp '+identifier+'.tgz $CMSSW_BASE/../')
    shell_finisher.close()