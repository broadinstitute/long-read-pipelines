TAG1 = kgarimella/gatk-lr-custom:0.01.17
TAG2 = kgarimella/gatk-lr-custom:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
