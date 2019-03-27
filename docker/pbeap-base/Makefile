TAG1 = kgarimella/pbeap-base:v9
TAG2 = kgarimella/pbeap-base:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
