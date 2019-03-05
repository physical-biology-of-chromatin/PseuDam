#/bin/sh

# A POSIX variable
OPTIND=1

# Initialize our own variables:
tool=""
version=""

while getopts "h?v:t:p:" opt; do
  case "$opt" in
  h|\?)
    echo "update_tools.sh -t toolname -v tool_version -p tool_previous_version"
    exit 0
    ;;
  v)
    version=$OPTARG
    ;;
  p)
    prev_version=$OPTARG
    ;;
  t)
    tool=$OPTARG
    ;;
  esac
done

echo "tool=$tool, version='$version'"

docker_tool_dir="src/docker_modules/"$tool"/"
echo $docker_tool_dir
if [ -d $docker_tool_dir ]; then
  echo "docker module found for $tool."
  if [ -d $docker_tool_dir$version ]; then
    echo "version already existing, skipping."
  else
    cp -R $docker_tool_dir$prev_version $docker_tool_dir$version
    sed -i "s|$prev_version|$version|g" "$docker_tool_dir$version/Dockerfile"
    sed -i "s|$prev_version|$version|g" "$docker_tool_dir$version/docker_init.sh"
    echo "docker_module for $tool:$version, done."
  fi
else
  echo "docker module not found for '$tool', skipping."
fi

singularity_tool_dir="src/singularity_modules/"$tool"/"
echo $singularity_tool_dir
if [ -d $singularity_tool_dir ]; then
  echo "singularity module found for $tool."
  if [ -d $singularity_tool_dir$version ]; then
    echo "version already existing, skipping."
  else
    cp -R $singularity_tool_dir$prev_version $singularity_tool_dir$version
    sed -i "s|$prev_version|$version|g" "$singularity_tool_dir$version/$tool.def"
    sed -i "s|$prev_version|$version|g" "$singularity_tool_dir$version/build.sh"
    echo "singularity_module for $tool:$version, done."
  fi
else
  echo "singularity module not found for '$tool', skipping."
fi

nf_tool_dir="src/nf_modules/"$tool"/"
echo $nf_tool_dir
if [ -d $nf_tool_dir ]; then
  echo "nf module found for $tool."
  find $nf_tool_dir -maxdepth 1 -mindepth 1 -type f -name "*.config" |
    awk "{system(\"sed -i \\\"s|$prev_version|$version|g\\\" \"\$0)}"
  echo "nf_module for $tool:$version, done."
else
  echo "nf module not found for '$tool', skipping."
fi
