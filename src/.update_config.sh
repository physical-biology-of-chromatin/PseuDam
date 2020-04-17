# update docker url
fd ".*config" -E "nf_modules" src/ -x perl -0777pe 's|container = "|container = "lbmc/|g' -i {}

# update singularity url
fd ".*config" -E "nf_modules" src/ -x perl -pe 's|container = "lbmc/file://bin/(.*).img"|container = "lbmc/\1"|g' -i {}

# update singularity config
fd ".*config" -E "nf_modules" src/ -x perl -0777pe 's|\n\s*singularity {\n\s*singularity.enabled = true|\n  singularity {\n    singularity.enabled = true\n    singularity.cacheDir = "./bin/"|mg' -i {}

# update in2p3 config
fd ".*config" -E "nf_modules" src/ -x perl -0777pe 's|\n\s*ccin2p3 {\n\s*singularity.enabled = true|\n  ccin2p3 {\n    singularity.enabled = true\n    singularity.cacheDir = "/sps/lbmc/common/singularity/"|mg' -i {}
fd ".*config" src/ -x perl -pe 's|container = "lbmc//sps/lbmc/common/singularity/(.*).img"|container = "lbmc/\1"|g' -i {}
fd ".*config" -E "nf_modules" src/ -x perl -0777pe 's|singularity.cacheDir = "/sps/lbmc/common/singularity/"|singularity.cacheDir = "\$baseDir/.singularity_in2p3/"|mg' -i {}

# we remove the ccin2p3_conda section
fd ".*config" -E "nf_modules" src/ -x perl -0777pe "s|\s*ccin2p3_conda {.*ccin2p3 {\n|\n  ccin2p3 {\n|msg" -i {}

# we update the psmn module to conda
fd ".*config" -E "nf_modules" src/ -x perl -0777pe 's|beforeScript = "source /usr/share/lmod/lmod/init/bash; module use ~/privatemodules"\n\s*module = "(.*)/(.*)"|beforeScript = "source \$baseDir/.conda_psmn.sh"\n        conda = "\$baseDir/.conda_envs/\L\1_\2"|mg' -i {}

# we update the psmn queue to new cluster
fd ".*config" src/ -x perl -0777pe 's|E5-2670deb128A,E5-2670deb128B,E5-2670deb128C,E5-2670deb128D,E5-2670deb128E,E5-2670deb128F|CLG6242deb384A,CLG6242deb384C,CLG5218deb192A,CLG5218deb192B,CLG5218deb192C,CLG5218deb192D,SLG5118deb96,SLG6142deb384A,SLG6142deb384B,SLG6142deb384C,SLG6142deb384D|mg' -i {}
fd ".*config" src/ -x perl -0777pe 's|monointeldeb128,monointeldeb48,h48-E5-2670deb128,h6-E5-2667v4deb128|monointeldeb128|mg' -i {}
fd ".*config" src/ -x perl -0777pe 's|openmp16|openmp32|mg' -i {}
fd ".*config" src/ -x perl -0777pe 's|cpus = 16|cpus = 32|mg' -i {}
