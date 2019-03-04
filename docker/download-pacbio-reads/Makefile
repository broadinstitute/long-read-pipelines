TAG1 = weisburd/download-pacbio-study:v1.4
TAG2 = weisburd/download-pacbio-study:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
