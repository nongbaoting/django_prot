# Generated by Django 2.2.5 on 2022-04-04 20:33

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0028_protcd_protcdncbi'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='protcdncbi',
            name='cdd_id',
        ),
        migrations.DeleteModel(
            name='protCD',
        ),
        migrations.DeleteModel(
            name='protCDncbi',
        ),
    ]