# Generated by Django 2.2.5 on 2022-04-04 18:03

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0025_nrinfo_protcd_protcdncbi'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='protcd',
            name='cdd_name',
        ),
        migrations.RemoveField(
            model_name='protcdncbi',
            name='cdd_name',
        ),
    ]