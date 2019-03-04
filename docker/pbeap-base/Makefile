TAG1 = kgarimella/pbeap-base:v1
TAG2 = kgarimella/pbeap-base:latest

all: build tag push

build:
	docker build -t $(TAG1) -t $(TAG2) .

tag:
	docker tag pbeap-base $(TAG1)
	docker tag pbeap-base $(TAG2)

push:
	docker push $(TAG1)
	docker push $(TAG2)
