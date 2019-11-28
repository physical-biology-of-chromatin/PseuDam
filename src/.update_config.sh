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
