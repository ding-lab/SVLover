export PATH=$(echo $PATH | tr ':' '\n' | grep -v '.sdkman' | tr '\n' ':')
JAVA_CMD: /home/yuweiz/.sdkman/candidates/java/current/bin/java
JAVA_HOME: /home/yuweiz/.sdkman/candidates/java/current
unset JAVA_CMD
unset JAVA_HOME

# 2. (Optional but recommended) Verify that your Conda Java is the default one again:
echo "Current Java Path being used by shell: $(which java)"

# 3. Re-run Nextflow:
# nextflow run main.nf -c nextflow.config

nextflow -log nextflow.log  run main.nf -c nextflow.config

# docker run -it wwliao/svision-pro:2.4--f19eb78 /bin/bash