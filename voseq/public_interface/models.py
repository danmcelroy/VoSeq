from django.db import models


class Genes(models.Model):
    gene_code = models.CharField(max_length=255)
    genetic_code = models.PositiveSmallIntegerField(blank=True)
    length = models.PositiveSmallIntegerField()
    description = models.CharField(max_length=255, blank=True)
    reading_frame = models.PositiveSmallIntegerField()
    notes = models.TextField(blank=True)
    aligned = models.CharField(
        max_length=6, choices=(
            ('yes', 'yes'),
            ('no', 'no'),
            ('notset', 'notset'),
        ),
        default='notset',
    )
    intron = models.CharField(max_length=255, blank=True)
    prot_code = models.CharField(
        max_length=6, choices=(
            ('yes', 'yes'),
            ('no', 'no'),
            ('notset', 'notset'),
        ),
        default='notset',
    )
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

    DONT_KNOW = 'd'
    YES = 'y'
    NO = 'n'
    TYPE_SPECIES_CHOICES = (
        (DONT_KNOW, 'don\'t know'),
        (YES, 'yes'),
        (NO, 'no'),
    )

    SPREAD = 's'
    ENVELOPE = 'e'
    PHOTO = 'p'
    NONE = 'n'
    DESTROYED = 'd'
    LOST = 'l'
    VOUCHER_CHOICES = (
        (SPREAD, 'spread'),
        (ENVELOPE, 'in envelope'),
        (PHOTO, 'only photo'),
        (NONE, 'no voucher'),
        (DESTROYED, 'destroyed'),
        (LOST, 'lost'),
    )
    code = models.CharField(max_length=100, unique=True,
                            primary_key=True,
                            help_text="Voucher code.",
                            )
    orden = models.CharField(max_length=100, blank=True, null=True)
    superfamily = models.CharField(max_length=100, null=True)
    family = models.CharField(max_length=100, blank=True, null=True)
    subfamily = models.CharField(max_length=100, blank=True, null=True)
    tribe = models.CharField(max_length=100, blank=True, null=True)
    subtribe = models.CharField(max_length=100, blank=True, null=True)
    genus = models.CharField(max_length=100, blank=True, null=True)
    species = models.CharField(max_length=100, blank=True, null=True)
    subspecies = models.CharField(max_length=100, blank=True)
    country = models.CharField(max_length=100, blank=True, null=True)
    specificLocality = models.CharField(max_length=250, help_text="Locality of origin for this specimen.", blank=True, null=True)
    typeSpecies = models.CharField(max_length=1, choices=TYPE_SPECIES_CHOICES,
                                   help_text="Is this a type species?", blank=True, null=True)
    latitude = models.FloatField(blank=True, null=True)
    longitude = models.FloatField(blank=True, null=True)
    max_altitude = models.IntegerField(blank=True, null=True, help_text="Enter altitude in meters above sea level.")
    min_altitude = models.IntegerField(blank=True, null=True, help_text="Enter altitude in meters above sea level.")
    collector = models.CharField(max_length=100, blank=True, null=True)
    dateCollection = models.DateField(null=True)  # TODO check if better blank null rather than null true
    extraction = models.CharField(max_length=50, help_text="Number of extraction event.", blank=True, null=True)
    extractionTube = models.CharField(max_length=50, help_text="Tube containing DNA extract.", blank=True, null=True)
    dateExtraction = models.DateField(null=True)
    extractor = models.CharField(max_length=100, blank=True, null=True)
    voucherLocality = models.CharField(max_length=200, blank=True, null=True)
    publishedIn = models.TextField(blank=True, null=True)
    notes = models.TextField(blank=True, null=True)
    edits = models.TextField(blank=True, null=True)
    latesteditor = models.TextField(blank=True, null=True)
    hostorg = models.CharField(max_length=200, help_text="Hostplant or other host.", blank=True, null=True)
    sex = models.CharField(max_length=1, choices=SEX_CHOICES, blank=True, null=True)
    voucher = models.CharField(max_length=1, choices=VOUCHER_CHOICES, blank=True, null=True)
    voucherCode = models.CharField(max_length=100, help_text="Original code of voucher specimen.", blank=True, null=True)
    determinedBy = models.CharField(max_length=100, help_text="Person that identified the taxon for this specimen.", blank=True, null=True)
    auctor = models.CharField(max_length=100, help_text="Person that described this taxon.", blank=True, null=True)
    timestamp = models.DateTimeField(null=True)  # TODO change this to date_created = models.DateField(auto_now_add=True)


class Sequences(models.Model):
    code = models.ForeignKey(Vouchers, help_text='Save as lower case.')
    gene_code = models.CharField(max_length=100)
    sequences = models.TextField(blank=True)
    accession = models.CharField(max_length=100, blank=True)
    labPerson = models.CharField(max_length=100, blank=True)
    time_created = models.DateTimeField(auto_now_add=True, null=True, blank=True)
    time_edited = models.DateTimeField(auto_now=True, null=True, blank=True)
    notes = models.TextField(blank=True)
    genbank = models.NullBooleanField()


class FlickrImages(models.Model):
    voucher = models.ForeignKey(
        Vouchers,
        help_text='Relation with id of voucher. Save as lower case.',
        blank=True,
        null=True,
    )
    voucherImage = models.URLField(help_text="URLs of the Flickr page.")
    thumbnail = models.URLField(help_text="URLs for the small sized image from Flickr.")
    flickr_id = models.CharField(max_length=100, help_text="ID numbers from Flickr for our photo.")
