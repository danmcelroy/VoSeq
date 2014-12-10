from django.db import models


class Genes(models.Model):
    gene_code = models.CharField(max_length=255)
    genetic_code = models.PositiveSmallIntegerField(blank=True)
    length = models.PositiveSmallIntegerField()
    description = models.CharField(max_length=255, blank=True)
    reading_frame = models.PositiveSmallIntegerField()
    notes = models.TextField(blank=True)
    aligned = models.CharField(max_length=6, choices=(
        ('yes', 'yes'),
        ('no', 'no'),
        ('notset', 'notset'),
    ),
                               default='notset')
    intron = models.CharField(max_length=255, blank=True)
    prot_code = models.CharField(max_length=6, choices=(
        ('yes', 'yes'),
        ('no', 'no'),
        ('notset', 'notset'),
    ),
                               default='notset')
    gene_type = models.CharField(max_length=255, blank=True)

    time_created = models.DateTimeField(auto_now_add=True)


class GeneSets(models.Model):
    geneset_name = models.CharField(max_length=75, default=None)
    geneset_creator = models.CharField(max_length=75, default=None)
    geneset_description = models.CharField(max_length=100, default=None, blank=True)
    geneset_list = models.TextField(default=None)


class Members(models.Model):
    firstname = models.CharField(max_length=100)
    lastname = models.CharField(max_length=100)
    login = models.CharField(max_length=100)
    passwd = models.CharField(max_length=100)
    admin = models.BinaryField(default=None)


class Primers(models.Model):
    code = models.CharField(max_length=100)
    gene_code = models.CharField(max_length=100)
    primer1 = models.CharField(max_length=100, blank=True)
    primer2 = models.CharField(max_length=100, blank=True)
    primer3 = models.CharField(max_length=100, blank=True)
    primer4 = models.CharField(max_length=100, blank=True)
    primer5 = models.CharField(max_length=100, blank=True)
    primer6 = models.CharField(max_length=100, blank=True)


class Sequences(models.Model):
    code = models.CharField(max_length=100)
    gene_code = models.CharField(max_length=100)
    sequences = models.TextField()
    accession = models.CharField(max_length=100, blank=True)
    labPerson = models.CharField(max_length=100, blank=True)
    time_created = models.DateTimeField(auto_now_add=True)
    time_edited = models.DateTimeField(auto_now=True)
    notes = models.TextField(blank=True)
    genbank = models.BooleanField(default=None, blank=True)


class TaxonSets(models.Model):
    taxonset_name = models.CharField(max_length=50)
    taxonset_creator = models.CharField(max_length=75)
    taxonset_description = models.CharField(max_length=100)
    taxonset_list = models.TextField()


class Vouchers(models.Model):
    MALE = 'm'
    FEMALE = 'f'
    LARVA = 'l'
    WORKER = 'w'
    QUEEN = 'q'
    SEX_CHOICES = (
        (MALE, 'male'),
        (FEMALE, 'female'),
        (LARVA, 'larva'),
        (WORKER, 'worker'),
        (QUEEN, 'queen'),
    )

    SPREAD = 's'
    ENVELOPE = 'e'
    PHOTO = 'p'
    NONE = 'n'
    VOUCHER_CHOICES = (
        (SPREAD, 'spread'),
        (ENVELOPE, 'in envelope'),
        (PHOTO, 'only photo'),
        (NONE, 'no voucher'),
    )
    code = models.CharField(max_length=100, unique=True, help_text="Voucher code.")
    orden = models.CharField(max_length=100)
    family = models.CharField(max_length=100)
    subfamily = models.CharField(max_length=100)
    tribe = models.CharField(max_length=100)
    subtribe = models.CharField(max_length=100)
    genus = models.CharField(max_length=100)
    species = models.CharField(max_length=100)
    subspecies = models.CharField(max_length=100)
    country = models.CharField(max_length=100)
    specificLocality = models.CharField(max_length=100, help_text="Locality of origin for this specimen.")
    typeSpecies = models.CharField(max_length=100, help_text="Is this a type species?")
    latitude = models.FloatField()
    longitude = models.FloatField()
    max_altitude = models.IntegerField(help_text="Enter altitude in meters above sea level.")
    min_altitude = models.IntegerField(help_text="Enter altitude in meters above sea level.")
    collector = models.CharField(max_length=100)
    dateCollection = models.DateField(null=True) # TODO check if better blank null rather than null true
    extraction = models.CharField(max_length=50, help_text="Number of extraction event.")
    extractionTube = models.CharField(max_length=50, help_text="Tube containing DNA extract.")
    dateExtraction = models.DateField(null=True)
    extractor = models.CharField(max_length=100)
    voucherLocality = models.CharField(max_length=200)
    publishedIn = models.TextField()
    notes = models.TextField()
    edits = models.TextField()
    latesteditor = models.TextField()
    hostorg = models.CharField(max_length=200, help_text="Hostplant or other host.")
    sex = models.CharField(max_length=1, choices=SEX_CHOICES)
    voucher = models.CharField(max_length=1, choices=VOUCHER_CHOICES)
    voucherCode = models.CharField(max_length=100, help_text="Original code of voucher specimen.")
    determinedBy = models.CharField(max_length=100, help_text="Person that identified the taxon for this specimen.")
    auctor = models.CharField(max_length=100, help_text="Person that described this taxon.")
    timestamp = models.DateTimeField()

class FlickrImages(models.Model):
    voucher = models.ForeignKey('Vouchers', help_text='Relation with id of voucher')
    voucherImage = models.URLField(help_text="URLs of the Flickr page.")
    thumbnail = models.URLField(help_text="URLs for the small sized image from Flickr.")
    flickr_id = models.IntegerField(help_text="ID numbers from Flickr for our photo.")
