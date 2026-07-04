IMAGE_NAME = lr-flare
FLARE_VERSION = 0.6.0
IMAGE_VERSION = v3
VERSION = $(FLARE_VERSION)-$(IMAGE_VERSION)

TAG1 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):$(VERSION)
TAG2 = us.gcr.io/broad-dsp-lrma/$(IMAGE_NAME):latest

all: | build push

build:
	docker build --platform linux/amd64 -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache --platform linux/amd64 -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)

# Build natively on amd64 in GCP and push to GCR (recommended on Apple Silicon).
cloudbuild:
	gcloud builds submit --project broad-dsp-lrma --config cloudbuild.yaml .
