# Generated by Django 2.2.5 on 2022-04-08 20:52

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0043_protcdncbi'),
    ]

    operations = [
        migrations.RenameField(
            model_name='submitinfonew',
            old_name='proj_name',
            new_name='job_name',
        ),
        migrations.RenameField(
            model_name='submitinfonew',
            old_name='run_status',
            new_name='result',
        ),
    ]
