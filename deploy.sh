#!/bin/bash
cd ~/project/
commitID=`git log -n 1 --pretty="%h" -- environment.yml`
sed -i '/^# environment.yml/d' Singularity && echo -e "\n# environment.yml commit ID: $commitID\n" >> Singularity
git config --global user.email "alcalan@iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add dag.png
git add dag.html
git status
git commit -m "circle CI deployment [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://hub.docker.com/api/build/v1/source/f121de2c-2362-489f-b8d6-d57fa8a23685/trigger/452409e9-44f3-4b82-959d-8058b40116c1/call/
