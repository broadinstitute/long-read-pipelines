FROM weisburd/base-image-for-str-tools:v1

MAINTAINER Ben Weisburd

COPY GATK_withSTRmodel.jar /
COPY wgs-all.tab /

# remove ENTRYPOINT because cromwell doesn't work with it defined
#ENTRYPOINT [ "java", "-jar", "GATK_withSTRmodel.jar" ]