# Generated by Django 2.2.5 on 2022-04-05 19:35

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0034_auto_20220405_0824'),
    ]

    operations = [
        migrations.CreateModel(
            name='protCDncbiOne',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('cdd_nameCat', models.TextField()),
                ('cdd_idCat', models.TextField()),
                ('protin_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='nrinfo_protCDncbiOne', to='protein.NrInfo')),
            ],
        ),
    ]