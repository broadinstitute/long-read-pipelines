IMAGE_NAME=lr-ontqc
VERSION=1.0.0

TAG1="us-central1-docker.pkg.dev/broad-dsp-lrma/general/$(IMAGE_NAME):$(VERSION)"
TAG2="us-central1-docker.pkg.dev/broad-dsp-lrma/general/$(IMAGE_NAME):latest"

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
